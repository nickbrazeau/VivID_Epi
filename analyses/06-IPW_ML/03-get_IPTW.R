#----------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate the IPTW
#----------------------------------------------------------------------------------------------------

# imports
library(tidyverse)
source("R/00-functions_iptw.R")




# read in param table and results of training
params <- readRDS("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/final_models/_rslurm_vivid_preds/params.RDS")
trainpaths <-  list.files(path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/final_models/_rslurm_vivid_preds/", pattern = ".RDS", full.names = T)
trainpaths <- trainpaths[!c(grepl("params.RDS", trainpaths) | grepl("f.RDS", trainpaths))]

# sort properly to match rows in df
trainpaths <- tibble::tibble(trainpaths = trainpaths) %>% 
  mutate(order = stringr::str_extract(basename(trainpaths), "[0-9]+"),
         order = as.numeric(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(-c(order)) %>% 
  unlist(.)
params$trained <- purrr::map(  trainpaths, function(x){ ret <- readRDS(x); return(ret[[1]]) }  )

# get predictions
params$preds <- purrr::pmap(params[,c("task", "trained")], function(task, trained){ ret <- predict(trained, task); return(ret) })


# make inverse probability weights
params$iptw <- pmap(params[,c("task", "preds")], get_iptw_prob)



#................................................
# Analyze Weights
#................................................
lapply(params$iptw, summary)












