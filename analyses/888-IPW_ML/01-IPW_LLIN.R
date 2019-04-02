#----------------------------------------------------------------------------------------------------
# Purpose of this to calculate IPW for LLIN
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(srvyr)
library(mlr)
source("R/00-functions_basic.R")
set.seed(42)

#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "internal_datamap_files/DERIVED_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model),
                llin_causal_model = ifelse(is.na(llin_causal_model), "n", llin_causal_model)
  )

sf::st_geometry(dt) <- NULL


#......................
# DAG
#......................
# Based on our DAG (http://dagitty.net/development/dags.html?id=jKNQhB#)
# we have identified wealth, education, age, and sex as confounders 
# but will use all potential confounders
llinmdl <- dcdr$column_name[dcdr$llin_causal_model == "y"]

llin <- dt %>% 
  dplyr::select(c("pv18s" , "pfldh", llinmdl))

#......................
# final preprocess
#......................
# check for class imbalance
xtabs(~llin$hml20_fctb)
# make sure recoding worked
mlr::summarizeColumns(llin) %>% 
  dplyr::mutate_if(is.numeric, round, 2)
# check missingness
mice::md.pattern(llin, rotate.names = T)

# MCAR assumption (revisit)
llin <- llin %>% 
  dplyr::filter(complete.cases(.))



#-------------------------------------------------------------------
##########                 Ensemble Setup                 ##########
#-------------------------------------------------------------------
# RESAMPLING 
resamp_sp <- mlr::makeResampleDesc("SpRepCV", fold = 5, reps = 100)

# TRAIN set & TEST set
n <- nrow(llin)
partition <- sample(1:n, size = n*0.8)
full <-   llin[,           !colnames(llin) %in% c("pfldh", "pv18s", "hv005_wi")]
train <-  llin[partition,  !colnames(llin) %in% c("pfldh", "pv18s", "hv005_wi")]
test <-   llin[-partition, !colnames(llin) %in% c("pfldh", "pv18s", "hv005_wi")]

# TASK
llin.classif.train <- mlr::makeClassifTask(id = "LLIN_IPW", 
                                           data = train, 
                                           target = "hml20_fctb",
                                           positive = "yes",
                                           coordinates = train[,c("longnum", "latnum")])

llin.classif.test <- mlr::makeClassifTask( id = "LLIN_IPW", 
                                           data = test, 
                                           target = "hml20_fctb",
                                           positive = "yes",
                                           coordinates = test[,c("longnum", "latnum")])

llin.classif <- mlr::makeClassifTask(      id = "LLIN_IPW", 
                                           data = full, 
                                           target = "hml20_fctb",
                                           positive = "yes",
                                           coordinates = llin[,c("longnum", "latnum")])


#-------------------------------------------------------------------
##########                 Make Ensebmble                 ##########
#-------------------------------------------------------------------
# https://mlr.mlr-org.com/articles/tutorial/integrated_learners.html
base.lrns <- lapply(c("classif.logreg",
                      "classif.glmnet", 
#                      "classif.svm",
#                      "classif.gausspr",
#                      "classif.gbm",
#                      "classif.kknn",
                      "classif.randomForest"
),
makeLearner, predict.type = "prob")

stckd.lrn = makeStackedLearner(base.learners = base.lrns,
                               predict.type = "prob", method = "hill.climb")

mlr::getLearnerParamSet(stckd.lrn)

#-------------------------------------------------------------------
##########                 Tune Ensebmble                 ##########
#-------------------------------------------------------------------

# 
# # tune S
# num_ps = makeParamSet(
#   makeNumericParam("s", lower = 0.01, upper = 0.02) #, trafo = function(x){10^x})
# )
# ctrl <- makeTuneControlGrid()
# 
# glmnet.tuned <- mlr::tuneParams(learner = glmnet.lrn,
#                                 task = llin.classif,
#                                 resampling = resamp_sp,
#                                 measures = list(acc, auc),
#                                 par.set = num_ps,
#                                 control = ctrl)
# 
# glmnet.lrn = setHyperPars(makeLearner("classif.glmnet",
#                                       predict.type = "prob"), 
#                           par.vals = glmnet.tuned$x)
# 
# 
# # evaluate
# glmnet.perform <- mlr::resample(learner = glmnet.lrn,
#                                 task = llin.classif, 
#                                 resampling = resamp_sp,
#                                 measures = list(acc, auc, mmce),
#                                 keep.pred = F)
# make full model
stckd.mod <- train(learner = stckd.lrn,
                    task = llin.classif)


#-------------------------------------------------------------------
##########            Ensebmble Performance                ##########
#-------------------------------------------------------------------
# # make predictions
# glmnet.pred <- predict(glmnet.mod, 
#                        task = llin.classif.test)
#                             
# # test performance
# mlr::calculateConfusionMatrix(glmnet.pred)
# mlr::calculateROCMeasures(glmnet.pred)


#-------------------------------------------------------------------
##########            Ensebmble Predictions                ##########
#-------------------------------------------------------------------
# predict on full model
stckd.full <- predict(stckd.mod, 
                      task = llin.classif)
llin <- llin %>% 
  dplyr::mutate(hml20_b = ifelse(hml20_fctb == "no", 0, 1), # conver to binary
                pexp = mean(hml20_b)
  ) %>% 
  dplyr::bind_cols(., stckd.full$data) %>%
  dplyr::mutate(iptw_u = ifelse(hml20_fctb == "yes",
                                1/prob.yes,
                                1/(1-prob.yes)),
                iptw_s = pexp * iptw_u,
                iptwipsw = iptw_s * hv005_wi
                )

#-------------------------------------------------------------------
##########                Check Weights                   ##########
#-------------------------------------------------------------------
# check weight stability
summary(llin$iptw_u)
summary(llin$iptw_s)

# check weight distribution
weighted_llin <- llin %>% srvyr::as_survey_design(weights = iptwipsw)

llintable_orig <- tableone::CreateTableOne(vars = llinmdl[!llinmdl %in% c("pv18s", "pfldh", "hml20_fctb", "hv005_wi", "latnum", "longnum")], 
                                           strata = "hml20_fctb", test = F, 
                                           data = llin)

llintable_iptw <- tableone::svyCreateTableOne(vars = llinmdl[!llinmdl %in% c("pv18s", "pfldh", "hml20_fctb", "hv005_wi", "latnum", "longnum")], 
                                              strata = "hml20_fctb", test = F, 
                                              data = weighted_llin)


#-------------------------------------------------------------------
##########                Get Effects                     ##########
#-------------------------------------------------------------------
pv.llin.cest <- llin %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pv18s ~ hml20_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "logit"))

# broom::tidy(pv.llin.cest, exponentiate = T, conf.int = T)

pf.llin.cest <- llin %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pfldh ~ hml20_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "logit"))

# broom::tidy(pf.llin.cest, exponentiate = T, conf.int = T)


pv.llin.cest.tidy <- broom::tidy(pv.llin.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. vivax")

pf.llin.cest.tidy <- broom::tidy(pf.llin.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. falciparum")

pan.llin.cest.tidy <- dplyr::bind_rows(pv.llin.cest.tidy, pf.llin.cest.tidy) %>% 
  dplyr::select(c("species", "term", "estimate", "conf.low", "conf.high")) %>% 
  dplyr::filter(term != "(Intercept)") %>% 
  magrittr::set_colnames(c("species", "term", "OR", "L95", "U95")) 



save(llin, 
     llintable_orig, llintable_iptw,
     pv.llin.cest, pf.llin.cest,
     pan.llin.cest.tidy,
     file = "results/llin_iptw_models.rda")





