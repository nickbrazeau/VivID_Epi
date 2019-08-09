#----------------------------------------------------------------------------------------------------
# Purpose of this script is to apply Tuning Parameters and Make Final Models 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)
source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")

# after figuring out who won
# this will get hyperpar mlr::getHyperPars(learner.hp)

#............................................
# Apply the Tuning Results to the Learner
#............................................
params <- readRDS("analyses/06-IPW_ML/_rslurm_vivid_tunes/params.RDS")
# now overwrite Learner ModelMultiplexer to a 
# STACKED learner that is actually of an ensemble 
# 
# 
# 
#  CHANGED THIS TO AVG STACK
params$learner <- purrr::map(params$task, make_avg_Stack, 
                             learners = baselearners.list)
tuneresultpaths <- list.files(path = "analyses/06-IPW_ML/_rslurm_vivid_tunes/", pattern = ".RDS", full.names = T)
tuneresultpaths <- tuneresultpaths[!c(grepl("params.RDS", tuneresultpaths) | grepl("f.RDS", tuneresultpaths))]

# sort properly to match rows in df
tuneresultpaths <- tibble::tibble(tuneresultpaths = tuneresultpaths) %>% 
  mutate(order = stringr::str_extract(basename(tuneresultpaths), "[0-9]+"),
         order = as.numeric(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(-c(order)) %>% 
  unlist(.)


params$tuneresult <- purrr::map(tuneresultpaths, findbesttuneresult.simple)
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
                            jobname = 'vivid_preds_finalmodels',
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
























