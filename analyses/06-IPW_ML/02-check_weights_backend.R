#----------------------------------------------------------------------------------------------------
# Purpose of this script is to check for 
# covariate balance after applying the weights
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(energy)
source("R/00-IPTW_functions.R")

#............................................................
# Import Data
#............................................................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")

params <- readRDS("analyses/06-IPW_ML/_rslurm_vivid_spSL/params.RDS")
SLpaths <- list.files("analyses/06-IPW_ML/_rslurm_vivid_spSL/", 
                      pattern = ".RDS", full.names = T)

SLpaths <- SLpaths[!c(grepl("params.RDS", SLpaths) | grepl("f.RDS", SLpaths))]

# sort properly to match rows in df
SLpaths <- tibble::tibble(SLpaths = SLpaths) %>% 
  mutate(order = stringr::str_extract(basename(SLpaths), "[0-9]+"),
         order = as.numeric(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(-c(order)) %>% 
  unlist(.)
params$ensembl_cvRisk <- purrr::map(  SLpaths, function(x){ ret <- readRDS(x); return(ret[[1]]) }  )

#............................................................
# Get IPTW Estimates
#............................................................
params$SLpreds <-purrr::map(params$ensembl_cvRisk, "SL.predictions")
params$iptw <- purrr::pmap(params[,c("task", "SLpreds")], get_iptw_prob)

#............................................................
# Make Model Map
#............................................................
select_samples_by_iptw_and_feattarget <- function(target, feat, iptw, data){
  wi <- sample(1:nrow(data), prob = iptw, replace = T)
  ret <- data %>% 
    dplyr::select(c(target, feat))
  
  ret <- ret[wi, ]
  return(ret)
  
}


#............................................................
# Note, I want to sample 100 realizations of my
# sample weights since I can't apply the weights directly
# in the energy function 
#............................................................
modelmap.ml <- params %>% 
  dplyr::mutate(target = purrr::map(.$task, mlr::getTaskTargetNames)) %>% 
  dplyr::select(c("task", "iptw", "target")) %>% 
  tidyr::unnest(col = target)

modelmap.ml <- replicate(n = 100, modelmap.ml,
                         simplify = F) %>% 
  dplyr::bind_rows()

modelmap.ml$feat <- map(modelmap.ml$task, mlr::getTaskFeatureNames)
modelmap.ml$data <- lapply(1:nrow(modelmap.ml), function(x) return(dt))

# unnest to get all feats
modelmap.ml <- modelmap.ml %>% 
  tidyr::unnest(cols = feat)



modelmap.ml$data <- purrr::pmap(modelmap.ml[,c("target", "feat", "iptw", "data")], 
                                select_samples_by_iptw_and_feattarget)

modelmap.ml <- modelmap.ml %>% 
  dplyr::select(c("target", "feat", "data")) %>% 
  dplyr::rename(covar1 = target,
                covar2 = feat) %>% 
  dplyr::mutate(lvl = "IPTW")


# get base target-feats too
modelmap.nowi <- params %>% 
  dplyr::mutate(target = purrr::map(.$task, mlr::getTaskTargetNames)) %>% 
  dplyr::select(c("task", "target")) %>% 
  tidyr::unnest(col = target)

modelmap.nowi$feat <- map(modelmap.nowi$task, mlr::getTaskFeatureNames)
modelmap.nowi$data <- lapply(1:nrow(modelmap.nowi), function(x) return(dt))

modelmap.nowi <- modelmap.nowi %>% 
  tidyr::unnest(col = feat) %>% 
  dplyr::select(c("target", "feat", "data")) %>% 
  dplyr::rename(covar1 = target,
                covar2 = feat) %>% 
  dplyr::mutate(lvl = "base")



modelmap <- rbind.data.frame(modelmap.nowi, modelmap.ml)


#............................................................
# Energy function for slurm
#............................................................

slurm_calc_corr <- function(covar1, covar2, data, lvl){
  x1 <- unlist( data[,covar1] )
  x2 <- unlist( data[,covar2] )
  if(is.factor(x1)){
    x1 <- as.numeric(x1)
  }
  if(is.factor(x2)){
    x2 <- as.numeric(x2)
  }
  
  ret <- energy::dcor(x = x1, y = x2)
  
  return(ret)
}


#................................
# send it out with rSLURM
#................................

# for slurm on LL
setwd("analyses/06-IPW_ML/")
sjob <- rslurm::slurm_apply(f = slurm_calc_corr, 
                            params = modelmap, 
                            jobname = 'covar_IPTW_balance',
                            nodes = 1028, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 1028, 
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5:00:00"))
