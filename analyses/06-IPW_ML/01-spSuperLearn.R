#----------------------------------------------------------------------------------------------------
# Purpose of this script is to Ob-train the 
# Cross Validated Risk 
#----------------------------------------------------------------------------------------------------
remotes::install_github("nickbrazeau/mlrwrapSL"); library(mlrwrapSL)
library(tidyverse)
library(mlr)
library(rslurm)
source("R/00-functions_basic.R")
source("R/00-IPTW_functions.R")
source("analyses/06-IPW_ML/00-import_learners.R")
set.seed(48, "L'Ecuyer")

#...............................................................................................
# Import Data & tx map
#...............................................................................................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA))) %>% 
  dplyr::rename(target = column_name)

dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")

#...............................................................................................
# Read In Spatial Partition & make spatial cross-validation set
#...............................................................................................
drcpart <- readRDS("data/derived_data/kmean_drcpartitioned.rds")
dt.ml.cc <- dplyr::left_join(dt, drcpart, by = "hv001")
spcrossvalset <- split(1:nrow(dt.ml.cc), factor(dt.ml.cc$kmeansk))


#........................
# Subset and Store Dataframes for Tasks for prediciton on full dataset
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
# Manual processing from EDA, 
# know some of these parameters need to be changed
#...............................................................................................
txs$learnerlib[txs$target == "precip_mean_cont_scale_clst"] <- list(list(mlr::makeLearner("regr.lm", predict.type = "response")))
txs$learnerlib[txs$target == "hiv03_fctb"] <- list(list(mlr::makeLearner("classif.logreg", predict.type = "prob")))
txs$learnerlib[txs$target == "hv21345_fctb"] <- list(list(mlr::makeLearner("classif.logreg", predict.type = "prob")))
txs$learnerlib[txs$target == "ITN_fctb"] <- list(list(mlr::makeLearner("classif.logreg", predict.type = "prob")))



#...............................................................................................
# SUPERLEARNER
# Obtain Prediction Matrix
# Minimize Cross-validated Risk 
#...............................................................................................
paramsdf <- txs[,c("learnerlib", "task")]
paramsdf$valset.list <- lapply(1:nrow(paramsdf), function(x) return(spcrossvalset))



# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = mlrwrapSL::SL_crossval_risk_pred, 
                            params = paramsdf, 
                            jobname = 'vivid_spSL_final',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 17),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5-00:00:00"))

cat("*************************** \n Submitted SL models \n *************************** ")
























