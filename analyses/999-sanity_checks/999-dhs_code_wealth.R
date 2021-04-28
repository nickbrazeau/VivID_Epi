#----------------------------------------------------------------------------------------------------
# Purpose of this script is to run quick sensitivity of LLIN coding of nets
#----------------------------------------------------------------------------------------------------
remotes::install_github("nickbrazeau/mlrwrapSL"); library(mlrwrapSL)
library(tidyverse)
library(mlr)
library(energy)
options(survey.lonely.psu="adjust")
source("R/00-functions_basic.R")
source("R/00-IPTW_functions.R")
source("analyses/06-IPW_ML/00-import_learners.R")
set.seed(48, "L'Ecuyer")

#...............................................................................................
# Import Data & tx map
#...............................................................................................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map_v3.xlsx", sheet = 1) %>% 
  dplyr::pull(c("column_name"))
sf::st_geometry(dt) <- NULL


dt <- dt  %>% 
  dplyr::select(c("pv18s", "pfldh", "po18s", dcdr, "hv270_fctb")) %>% 
  dplyr::filter(complete.cases(.)) 


#............................................................
# base calc
#...........................................................
dtsrvy_raw <- dt %>% srvyr::as_survey_design(ids = hv001, 
                                             strata = hv023, 
                                             weights = hiv05_wi)
ret <- survey::svyglm("pv18s ~ hv270_fctb",
                      design = dtsrvy_raw,
                      family = quasibinomial(link="logit"))
broom::tidy(ret, exponentiate = T, conf.int = T)


ret <- survey::svyglm("pfldh ~ hv270_fctb",
                      design = dtsrvy_raw,
                      family = quasibinomial(link="logit"))
broom::tidy(ret, exponentiate = T, conf.int = T)


#............................................................
# Apply IPTW
#...........................................................
# read in treatments and subset to nets
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA))) %>% 
  dplyr::rename(target = column_name)
txs <- txs %>% 
  dplyr::filter(target == "wlthrcde_fctb") %>% 
  dplyr::mutate(target = "hv270_fctb",
                var_label = "DHS Wealth",
                abridged_var_label = "DHS Wealth",
                type = "binary")

#...............................................................................................
# Read In Spatial Partition & make spatial cross-validation set
#...............................................................................................
drcpart <- readRDS("data/derived_data/kmean_drcpartitioned.rds")
dt.ml.cc <- dplyr::left_join(dt, drcpart, by = "hv001")
spcrossvalset <- split(1:nrow(dt.ml.cc), factor(dt.ml.cc$kmeansk))


#........................
# Subset and Store Dataframes for Tasks for prediction on full dataset
#........................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml.cc %>% 
                            dplyr::select(c(x, y))
                          return(ret)
                        })


txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type")], 
                        .f = make_class_task)
#........................
# Make tasks and learner libraries
#........................
txs$learnerlib <- purrr::map(txs$type, function(x){
  if(x == "continuous"){
    return(base.learners.regr)
  } else if (x == "binary"){
    return(base.learners.classif)
  }
})


#...............................................................................................
# SUPERLEARNER
# Obtain Prediction Matrix
# Minimize Cross-validated Risk 
#...............................................................................................
paramsdf <- txs[,c("learnerlib", "task")]
paramsdf$valset.list <- lapply(1:nrow(paramsdf), function(x) return(spcrossvalset))

#......................
# run
#......................
paramsdf$ensembl_cvRisk <- purrr::pmap(paramsdf, mlrwrapSL::SL_crossval_risk_pred, 
                                       seed = 123)


#............................................................
# Look at distribution and effects of iptweights
#...........................................................
# tidy up
paramsdf <- paramsdf %>% 
  dplyr::mutate(target = purrr::map_chr(task, mlr::getTaskTargetNames)) %>% 
  dplyr::select(c("target", "ensembl_cvRisk"))

#......................
# Get IPTW Estimates
#......................
txs <- dplyr::left_join(txs, paramsdf, by = "target")
txs$SLpreds <- purrr::map(txs$ensembl_cvRisk, "SL.predictions")
txs$iptw <- purrr::pmap(txs[,c("task", "SLpreds")], get_iptw_prob)

# look at balance
boxplot(txs$iptw[[1]])
boxplot(log(txs$iptw[[1]]))
summary(txs$iptw[[1]])

#......................
# run again but now with simpler since
# urbanicity highly correlated and causing very high weights
#......................
txs <- txs %>% 
  dplyr::select(-c("ensembl_cvRisk", "SLpreds", "iptw"))
txs$learnerlib[txs$target == "hv270_fctb"] <- list(list(mlr::makeLearner("classif.logreg", predict.type = "prob")))
paramsdf <- txs[,c("learnerlib", "task")]
paramsdf$valset.list <- lapply(1:nrow(paramsdf), function(x) return(spcrossvalset))

# run
paramsdf$ensembl_cvRisk <- purrr::pmap(paramsdf, mlrwrapSL::SL_crossval_risk_pred, 
                                       seed = 123)


#............................................................
# Look at distribution and effects of iptweights
#...........................................................
# tidy up
paramsdf <- paramsdf %>% 
  dplyr::mutate(target = purrr::map_chr(task, mlr::getTaskTargetNames)) %>% 
  dplyr::select(c("target", "ensembl_cvRisk"))

#......................
# Get IPTW Estimates
#......................
txs <- dplyr::left_join(txs, paramsdf, by = "target")
txs$SLpreds <- purrr::map(txs$ensembl_cvRisk, "SL.predictions")
txs$iptw <- purrr::pmap(txs[,c("task", "SLpreds")], get_iptw_prob)

# look at balance
summary(txs$iptw[[1]])
boxplot(txs$iptw[[1]])
boxplot(log(txs$iptw[[1]]))



#............................................................
# apply in mod
#...........................................................
dt$wi <- unlist(txs$iptw)*dt$hiv05_wi
options(survey.lonely.psu="adjust")
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, 
                                         strata = hv023, 
                                         weights = wi)
ret <- survey::svyglm("pv18s ~ hv270_fctb",
                      design = dtsrvy,
                      family = quasibinomial(link="logit"))
broom::tidy(ret, exponentiate = T, conf.int = T)

ret <- survey::svyglm("pfldh ~ hv270_fctb",
                      design = dtsrvy,
                      family = quasibinomial(link="logit"))
broom::tidy(ret, exponentiate = T, conf.int = T)
