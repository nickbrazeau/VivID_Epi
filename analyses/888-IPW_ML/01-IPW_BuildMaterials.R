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
bldmtrlmdl <- dcdr$column_name[dcdr$bldmtrl_causal_model == "y"]

bldmtrl <- dt %>% 
  dplyr::select(c("pv18s" , "pfldh", "hv005_wi", bldmtrlmdl, 
                  "longnum", "latnum"))

#......................
# final preprocess
#......................
# check for class imbalance
xtabs(~bldmtrl$hv21345_fctb)
# make sure recoding worked
mlr::summarizeColumns(bldmtrl) %>% 
  dplyr::mutate_if(is.numeric, round, 2)
# check missingness
mice::md.pattern(bldmtrl, rotate.names = T)

# MCAR assumption (revisit)
bldmtrl <- bldmtrl %>% 
  dplyr::filter(complete.cases(.))



#-------------------------------------------------------------------
##########                 Ensemble Setup                 ##########
#-------------------------------------------------------------------
# RESAMPLING 
resamp_sp <- mlr::makeResampleDesc("SpRepCV", fold = 5, reps = 100)

# TRAIN set & TEST set
n <- nrow(bldmtrl)
partition <- sample(1:n, size = n*0.8)
full <-   bldmtrl[,           !colnames(bldmtrl) %in% c("pfldh", "pv18s", "hv005_wi")]
train <-  bldmtrl[partition,  !colnames(bldmtrl) %in% c("pfldh", "pv18s", "hv005_wi")]
test <-   bldmtrl[-partition, !colnames(bldmtrl) %in% c("pfldh", "pv18s", "hv005_wi")]

# TASK
bldmtrl.classif.train <- mlr::makeClassifTask(id = "BuildMaterials_IPW", 
                                              data = train, 
                                              target = "hv21345_fctb",
                                              positive = "modern",
                                              coordinates = train[,c("longnum", "latnum")])

bldmtrl.classif.test <- mlr::makeClassifTask( id = "BuildMaterials_IPW", 
                                              data = test, 
                                              target = "hv21345_fctb",
                                              positive = "modern",
                                              coordinates = test[,c("longnum", "latnum")])

bldmtrl.classif <- mlr::makeClassifTask(      id = "BuildMaterials_IPW", 
                                              data = full, 
                                              target = "hv21345_fctb",
                                              positive = "modern",
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
ovrsmpl <- sum(bldmtrl$hv21345_fctb == "traditional")/sum(bldmtrl$hv21345_fctb == "modern")
stckd.lrn.smote <- mlr::makeSMOTEWrapper(stckd.lrn, sw.rate = floor(ovrsmpl), sw.nn = 5)


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
                   task = bldmtrl.classif.train)


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
                      task = bldmtrl.classif)

bldmtrl <- bldmtrl %>% 
  dplyr::mutate(hv21345_b = ifelse(hv21345_fctb == "traditional", 0, 1), # conver to binary
                pexp = mean(hv21345_b)
  ) %>% 
  dplyr::bind_cols(., stckd.full$data) %>%
  dplyr::mutate(iptw_u = ifelse(hv21345_fctb == "modern",
                                1/prob.modern,
                                1/(1-prob.modern)),
                iptw_s = pexp * iptw_u,
                iptwipsw = iptw_s * hv005_wi
  )

#-------------------------------------------------------------------
##########                Check Weights                   ##########
#-------------------------------------------------------------------
# check weight stability
summary(bldmtrl$iptw_u)
summary(bldmtrl$iptw_s)

# check weight distribution
weighted_bldmtrl <- bldmtrl %>% srvyr::as_survey_design(weights = iptwipsw)

bldmtrltable_orig <- tableone::CreateTableOne(vars = bldmtrlmdl[!bldmtrlmdl %in% c("pv18s", "pfldh", "hv21345_fctb", "hv005_wi", "latnum", "longnum")], 
                                           strata = "hv21345_fctb", test = F, 
                                           data = bldmtrl)

bldmtrltable_iptw <- tableone::svyCreateTableOne(vars = bldmtrlmdl[!bldmtrlmdl %in% c("pv18s", "pfldh", "hv21345_fctb", "hv005_wi", "latnum", "longnum")], 
                                              strata = "hv21345_fctb", test = F, 
                                              data = weighted_bldmtrl)


#-------------------------------------------------------------------
##########                Get Effects                     ##########
#-------------------------------------------------------------------
pv.bldmtrl.cest <- bldmtrl %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pv18s ~ hv21345_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "log"))

# broom::tidy(pv.bldmtrl.cest, exponentiate = T, conf.int = T)

pf.bldmtrl.cest <- bldmtrl %>% 
  dplyr::mutate(id = 1:nrow(.)) %>% 
  geepack::geeglm(pfldh ~ hv21345_fctb,
                  weights = iptwipsw,
                  data = .,
                  id = id,
                  family = binomial(link = "log"))

# broom::tidy(pf.bldmtrl.cest, exponentiate = T, conf.int = T)


pv.bldmtrl.cest.tidy <- broom::tidy(pv.bldmtrl.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. vivax")

pf.bldmtrl.cest.tidy <- broom::tidy(pf.bldmtrl.cest, exponentiate = T, conf.int = T) %>% 
  dplyr::mutate(species = "P. falciparum")

pan.bldmtrl.cest.tidy <- dplyr::bind_rows(pv.bldmtrl.cest.tidy, pf.bldmtrl.cest.tidy) %>% 
  dplyr::select(c("species", "term", "estimate", "conf.low", "conf.high")) %>% 
  dplyr::filter(term != "(Intercept)") %>% 
  magrittr::set_colnames(c("species", "term", "OR", "L95", "U95")) 



save(bldmtrl, 
     bldmtrltable_orig, bldmtrltable_iptw,
     pv.bldmtrl.cest, pf.bldmtrl.cest,
     pan.bldmtrl.cest.tidy,
     file = "results/bldmtrl_iptw_models.rda")





