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
# Manipulate tx map
#...............................................................................................
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
varstoinclude <- c("pv18s" , "pfldh", "hv001", "hv005_wi", txs$target,
                   "hab1_cont_scale", "hv104_fctb", # need to add in covariates that don't have confounding ancestors but are needed elsewhere
                   "urban_rura_fctb", "alt_dem_cont_scale_clst",
                   "hiv03_fctb", # no longer considered risk factor bc too few observations
                   "longnum", "latnum")

#...............................................................................................
# Import Data
#...............................................................................................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL

# subset to complete cases
dt.ml.cc <- dt  %>% 
  dplyr::select(varstoinclude) %>% 
  dplyr::filter(complete.cases(.)) 


#...............................................................................................
# Read In Spatial Partition & make spatial cross-validation set
#...............................................................................................
drcpart <- readRDS("data/derived_data/kmean_drcpartitioned.rds")
dt.ml.cc <- dplyr::left_join(dt.ml.cc, drcpart, by = "hv001")
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
# SUPERLEARNER
# Obtain Prediction Matrix
# Minimize Cross-validated Risk 
#...............................................................................................
paramsdf <- txs[,c("learnerlib", "task")]
paramsdf$valset.list <- lapply(1:nrow(paramsdf), function(x) return(spcrossvalset))


# start <- Sys.time()
# ret <- furrr::future_pmap(paramsdf, mlrwrapSL::SL_crossval_risk_pred)
# end <- Sys.time()
# end - start

# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = , 
                            params = paramsdf, 
                            jobname = 'vivid_spSL',
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

cat("*************************** \n Submitted SL models \n *************************** ")



pmap(paramsdf[3,], mlrwrapSL::SL_crossval_risk_pred)





















