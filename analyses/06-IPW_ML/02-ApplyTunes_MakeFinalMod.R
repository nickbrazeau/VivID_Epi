#----------------------------------------------------------------------------------------------------
# Purpose of this script is to Obtrain the 
# Cross Validated Risk 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)
library(rslurm)
source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-IPTW_functions.R")
source("R/00-Ensemble_CrossValidRisk.R")
# this will get hyperpar mlr::getHyperPars(learner.hp)


#...............................................................................................
# Set up data matrix to predict on
#...............................................................................................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL

#........................
# manipulate tx map
#........................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA))) %>% 
  dplyr::rename(target = column_name)



#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
varstoinclude <- c("pv18s" , "pfldh", "hv005_wi", txs$target,
                   "alt_dem_cont_scale_clst", "hab1_cont_scale", "hv104_fctb", "wtrdist_cont_log_scale_clst", # need to add in covariates that don't have confounding ancestors but are needed elsewhere
                   "hiv03_fctb", # no longer considered risk factor bc too few observations
                   "longnum", "latnum")

dt.ml <- dt %>% 
  dplyr::select(varstoinclude)

# subset to complete cases
dt.ml.cc <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum")) %>% 
  data.frame(.)


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
base.learners.classif <- lapply(baselearners.list$classif, function(x) return(mlr::makeLearner(x, predict.type = "prob")))
base.learners.regr <- lapply(baselearners.list$regress, function(x) return(mlr::makeLearner(x, predict.type = "response")))
txs$learnerlib <- purrr::map(txs$type, function(x){
  if(x == "continuous"){
    return(base.learners.regr)
  } else if (x == "binary"){
    return(base.learners.classif)
  }
})



#...............................................................................................
# Pull and Apply Tuning Results
#...............................................................................................
tuneresultpaths <- list.files(path = "analyses/06-IPW_ML/_rslurm_vivid_tunes_train", pattern = ".RDS", full.names = T)
tuneresultpaths <- tuneresultpaths[!c(grepl("params.RDS", tuneresultpaths) | grepl("f.RDS", tuneresultpaths))]

# sort properly to match rows in df
tuneresultpaths <- tibble::tibble(tuneresultpaths = tuneresultpaths) %>%
  mutate(order = stringr::str_extract(basename(tuneresultpaths), "[0-9]+"),
         order = as.numeric(order)) %>%
  dplyr::arrange(order) %>%
  dplyr::select(-c(order))

tuneresultpaths$tuneresult <- purrr::map(tuneresultpaths$tuneresultpaths, findbesttuneresult)
# apply
txs$learnerlib <- purrr::map2(txs$learnerlib, tuneresultpaths$tuneresult,  tune_learner_library)

#...............................................................................................
# Train Each Algorithm
# Obtain Prediction Matrix
# Minimize Cross-validated Risk 
#...............................................................................................


txs$learnerlib <- purrr::map(txs$type, function(x){
  if(x == "continuous"){
    return(base.learners.regr)
  } else if (x == "binary"){
    return(base.learners.classif)
  }
})


txs$proptrainset <- 0.5
paramsdf <- txs[,c("learnerlib", "task", "proptrainset")]

slurm_function <- function(learnerlib, task, proptrainset){
  source("R/00-Ensemble_CrossValidRisk.R") # need to add get_preds function to environment
  return( ensemble_crossval_risk_pred(learnerlib, task, proptrainset) )
}

# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = slurm_function, 
                            params = paramsdf, 
                            jobname = 'vivid_EL',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 17),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))

cat("*************************** \n Submitted EL models \n *************************** ")

























