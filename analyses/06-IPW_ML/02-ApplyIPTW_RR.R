#----------------------------------------------------------------------------------------------------
# Purpose of this script is to apply IPTW weights
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)

#......................................
# Read in param table sent to slurm
#......................................
resultsdir <- "~/Documents/MountPoints/mountedMeshnick/Sandbox/VivID_Epi/_rslurm_vivid_preds/"

results.list <- list.files(resultsdir, pattern = "results", full.names = T)
results <- lapply(results.list, function(x){
  
  ret <- readRDS(x)
  names(ret) <- basename(x)
  return(ret)
}) %>% 
  dplyr::bind_rows(., .id = "retname")


params <- readRDS(file = paste0(resultsdir, "params.RDS"))



