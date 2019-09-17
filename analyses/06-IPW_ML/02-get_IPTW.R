#----------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate the IPTW
#----------------------------------------------------------------------------------------------------

# imports
library(mlr)
library(tidyverse)
library(ggridges)
source("R/00-IPTW_functions.R")

# read in param table and results of training
#params <- readRDS("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/final_models/_rslurm_vivid_preds/params.RDS")
#trainpaths <-  list.files(path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/final_models/_rslurm_vivid_preds/", pattern = ".RDS", full.names = T)

params <- readRDS("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/_rslurm_vivid_spSL/params.RDS")
ELpaths <- list.files("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/_rslurm_vivid_spSL/", 
                      pattern = ".RDS", full.names = T)

ELpaths <- ELpaths[!c(grepl("params.RDS", ELpaths) | grepl("f.RDS", ELpaths))]

# sort properly to match rows in df
ELpaths <- tibble::tibble(ELpaths = ELpaths) %>% 
  mutate(order = stringr::str_extract(basename(ELpaths), "[0-9]+"),
         order = as.numeric(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(-c(order)) %>% 
  unlist(.)
params$ensembl_cvRisk <- purrr::map(  ELpaths, function(x){ ret <- readRDS(x); return(ret[[1]]) }  )

#...............................................................................................
# Get IPTW Estimates
#...............................................................................................
params$ELpreds <-purrr::map(params$ensembl_cvRisk, "EL.predictions")
params$iptw <- purrr::pmap(params[,c("task", "ELpreds")], get_iptw_prob)


#................................................
# Analyze Weights
#................................................
widist <- purrr::map(params$iptw, summary) %>% 
  do.call("rbind.data.frame", .) %>% 
  dplyr::mutate_if(is.numeric, round, 2) %>% 
  magrittr::set_colnames(c("min", "1stquart", "median", "mean", "3rdqart", "max"))

params$target <- unlist( purrr::map(params$task, mlr::getTaskTargetNames) )

widist <- cbind.data.frame(params$target, widist)

# make plot
plotdf <- params %>% 
  dplyr::select(c("target", "iptw")) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(iptw < 15)

ggplot(plotdf, aes(x = iptw, y = target)) +
  geom_density_ridges() + theme_ridges() +
  scale_y_discrete(expand = c(0.01, 0)) + 
  scale_x_continuous(expand = c(0.01, 0)) 


plotdf %>% 
  dplyr::group_by(target) %>% 
  dplyr::summarise(min = min(iptw),
                   mean = mean(iptw),
                   max = max(iptw))






