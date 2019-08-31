#----------------------------------------------------------------------------------------------------
# Purpose of this script is to apply Tuning Parameters and Make Final Models 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)
source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")

# after figuring out who won
# this will get hyperpar mlr::getHyperPars(learner.hp)



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

dt.ml.coords <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
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

#........................
# Pull Out Coords
#........................
txs$coordinates <- lapply(1:nrow(txs), function(x) return(dt.ml.coords))


#--------------------------------------
# Setup tasks & base learners
#--------------------------------------
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type", "coordinates")], 
                        .f = make_class_task)


#............................................
# Apply the Tuning Results to the Learner &
# update the "training" data to be the full data task
#............................................
params <- readRDS("analyses/06-IPW_ML/_rslurm_vivid_tunes_train/params.RDS")
# overwrite params tasks to the new tasks that we want to train
# and predict the data on the full data set now
params$task <- txs$task


params$learner <- purrr::map(params$task, make_avg_Stack, 
                             learners = baselearners.list)
tuneresultpaths <- list.files(path = "analyses/06-IPW_ML/_rslurm_vivid_tunes_train/", pattern = ".RDS", full.names = T)
tuneresultpaths <- tuneresultpaths[!c(grepl("params.RDS", tuneresultpaths) | grepl("f.RDS", tuneresultpaths))]

# sort properly to match rows in df
tuneresultpaths <- tibble::tibble(tuneresultpaths = tuneresultpaths) %>% 
  mutate(order = stringr::str_extract(basename(tuneresultpaths), "[0-9]+"),
         order = as.numeric(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(-c(order)) %>% 
  unlist(.)


params$tuneresult <- purrr::map(tuneresultpaths, findbesttuneresult)
params$tunedlearner <- purrr::pmap(params[, c("learner", "task", "tuneresult")], 
                                   tune_stacked_learner)
                                   

 #........................................................................
# Train on all the Data with Stacked Learner 
#........................................................................

slurm_traindata <- function(tunedlearner, task){
  ret <- mlr::train(learner = tunedlearner, 
                    task = task)
  return(ret)
  
}

paramsdf <- params[, c("tunedlearner", "task")]

# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = slurm_traindata, 
                            params = paramsdf, 
                            jobname = 'vivid_preds_finalmodels_average',
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

cat("*************************** \n Submitted final models \n *************************** ")



sjobs <- pmap(paramsdf, slurm_traindata)






















