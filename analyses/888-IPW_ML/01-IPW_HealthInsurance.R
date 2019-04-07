#----------------------------------------------------------------------------------------------------
# Purpose of this to calculate IPW for Health Insurance
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
dcdr <- readxl::read_excel(path = "internal_datamap_files/DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model),
                llin_causal_model = ifelse(is.na(llin_causal_model), "n", llin_causal_model),
                hlthinsrnc_causal_model = ifelse(is.na(hlthinsrnc_causal_model), "n", hlthinsrnc_causal_model),
                bldmtrl_causal_model = ifelse(is.na(bldmtrl_causal_model), "n", bldmtrl_causal_model)
                
  )
sf::st_geometry(dt) <- NULL


#......................
# DAG
#......................
# Based on our DAG (http://dagitty.net/development/dags.html?id=jKNQhB#)
# we have identified wealth, education, age, and sex as confounders 
# but will use all potential confounders
hlthinsnmdl <- dcdr$column_name[dcdr$hlthinsrnc_causal_model == "y"]

hlthins <- dt %>% 
  dplyr::select(c("pv18s" , "pfldh", "hv005_wi", hlthinsnmdl, 
                  "longnum", "latnum"))

#......................
# final preprocess
#......................
# check for class imbalance
xtabs(~hlthins$hab481_fctb)
# make sure recoding worked
mlr::summarizeColumns(hlthins) %>% 
  dplyr::mutate_if(is.numeric, round, 2)
# check missingness
mice::md.pattern(hlthins, rotate.names = T)

# MCAR assumption (revisit)
hlthins <- hlthins %>% 
  dplyr::filter(complete.cases(.))



#-------------------------------------------------------------------
##########                 Ensemble Setup                 ##########
#-------------------------------------------------------------------
# RESAMPLING 
resamp_sp <- mlr::makeResampleDesc("SpRepCV", fold = 5, reps = 100)

# TRAIN set & TEST set
n <- nrow(hlthins)
partition <- sample(1:n, size = n*0.8)
full <-   hlthins[,           !colnames(hlthins) %in% c("pfldh", "pv18s", "hv005_wi")]
train <-  hlthins[partition,  !colnames(hlthins) %in% c("pfldh", "pv18s", "hv005_wi")]
test <-   hlthins[-partition, !colnames(hlthins) %in% c("pfldh", "pv18s", "hv005_wi")]

# TASK
hlthins.classif.train <- mlr::makeClassifTask(id = "HealthInsurance_IPW", 
                                              data = train, 
                                              target = "hab481_fctb",
                                              positive = "yes",
                                              coordinates = train[,c("longnum", "latnum")])

hlthins.classif.test <- mlr::makeClassifTask( id = "HealthInsurance_IPW", 
                                              data = test, 
                                              target = "hab481_fctb",
                                              positive = "yes",
                                              coordinates = test[,c("longnum", "latnum")])

hlthins.classif <- mlr::makeClassifTask(      id = "HealthInsurance_IPW", 
                                              data = full, 
                                              target = "hab481_fctb",
                                              positive = "yes",
                                              coordinates = full[,c("longnum", "latnum")])


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

#-------------------------------------------------------------------
##########           Deal with Class Imbalance            ##########
#-------------------------------------------------------------------
# https://mlr.mlr-org.com/articles/tutorial/over_and_undersampling.html#tuning-the-probability-threshold
ovrsmpl <- sum(hlthins$hab481_fctb == "no")/sum(hlthins$hab481_fctb == "yes")
stckd.lrn.smote <- mlr::makeSMOTEWrapper(stckd.lrn, 
                                         sw.rate = floor(ovrsmpl), 
                                         sw.nn = 5)


#-------------------------------------------------------------------
##########                 Tune Ensebmble                 ##########
#-------------------------------------------------------------------
mlr::getLearnerParamSet(stckd.lrn)

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
                   task = hlthins.classif.train)

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
                      task = hlthins.classif)
hlthins <- hlthins %>% 
  dplyr::mutate(hab481_b = ifelse(hab481_fctb == "no", 0, 1), # conver to binary
                pexp = mean(hab481_b)
  ) %>% 
  dplyr::bind_cols(., stckd.full$data) %>%
  dplyr::mutate(iptw_u = ifelse(hab481_fctb == "yes",
                                1/prob.yes,
                                1/(1-prob.yes)),
                iptw_s = pexp * iptw_u,
                iptwipsw = iptw_s * hv005_wi
  )

#-------------------------------------------------------------------
##########                Check Weights                   ##########
#-------------------------------------------------------------------
# check weight stability
summary(hlthins$iptw_u)
summary(hlthins$iptw_s)

# check weight distribution
weighted_hlthins <- hlthins %>% srvyr::as_survey_design(weights = iptwipsw)

hlthinstable_orig <- tableone::CreateTableOne(vars = hlthinsnmdl[!hlthinsnmdl %in% c("pv18s", "pfldh", "hab481_fctb", "hv005_wi", "latnum", "longnum")], 
                                           strata = "hab481_fctb", test = F, 
                                           data = hlthins)

hlthinstable_iptw <- tableone::svyCreateTableOne(vars = hlthinsnmdl[!hlthinsnmdl %in% c("pv18s", "pfldh", "hab481_fctb", "hv005_wi", "latnum", "longnum")], 
                                              strata = "hab481_fctb", test = F, 
                                              data = weighted_hlthins)


#-------------------------------------------------------------------
##########                Get Effects                     ##########
#-------------------------------------------------------------------
pv.hlthins.cest <- hlthins %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pv18s ~ hab481_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "log"))

# broom::tidy(pv.hlthins.cest, exponentiate = T, conf.int = T)

pf.hlthins.cest <- hlthins %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pfldh ~ hab481_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "log"))

# broom::tidy(pf.hlthins.cest, exponentiate = T, conf.int = T)


pv.hlthins.cest.tidy <- broom::tidy(pv.hlthins.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. vivax")

pf.hlthins.cest.tidy <- broom::tidy(pf.hlthins.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. falciparum")

pan.hlthins.cest.tidy <- dplyr::bind_rows(pv.hlthins.cest.tidy, pf.hlthins.cest.tidy) %>% 
  dplyr::select(c("species", "term", "estimate", "conf.low", "conf.high")) %>% 
  dplyr::filter(term != "(Intercept)") %>% 
  magrittr::set_colnames(c("species", "term", "OR", "L95", "U95")) 



save(hlthins, 
     hlthinstable_orig, hlthinstable_iptw,
     pv.hlthins.cest, pf.hlthins.cest,
     pan.hlthins.cest.tidy,
     file = "results/hlthins_iptw_models.rda")





