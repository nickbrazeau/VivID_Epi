#----------------------------------------------------------------------------------------------------
# Purpose of this script is to apply Tuning Parameters and Make Final Models 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)
source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")

#......................................
# Read in param table sent to slurm
#......................................
resultsdir <- "analyses/06-IPW_ML/tune_modelparams/_rslurm_vivid_preds/"

results.list <- list.files(resultsdir, pattern = "results", full.names = T)
# need to do a proper sort
results.df <- tibble::tibble(name = basename(results.list),
                         path = results.list) %>% 
  dplyr::mutate(name = stringr::str_extract(name, "[0-9]+"),
                name = as.numeric(name)) %>% 
  dplyr::arrange(name)


results <- lapply(results.df$path, function(x){
  ret <- readRDS(x)
  return(ret[[1]])
}) 

params <- readRDS(file = paste0(resultsdir, "params.RDS"))
params$tuneresult <- lapply(1:length(results), function(x){
  return(results[[x]])
  })

#............................................
# Apply the Tuning Results to the Learner
#............................................
params$tunedlearner <- purrr::map(params$task, make_hillclimb_Stack, 
                          learners = baselearners.list)
params$tunedlearner <- purrr::pmap(params[, c("learner", "task", "tuneresult")], tune_stacked_learner)
                                   

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
setwd("analyses/06-IPW_ML/final_models")
ntry <- 18
sjob <- rslurm::slurm_apply(f = slurm_traindata, 
                            params = paramsdf, 
                            jobname = 'vivid_preds',
                            nodes = 18, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 128000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry - 1, 
                                                                 16),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))


cat("*************************** \n Submitted final models \n *************************** ")
























