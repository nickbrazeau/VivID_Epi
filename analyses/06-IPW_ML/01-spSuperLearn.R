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
source("analyses/06-IPW_ML/00-import_learners_updated.R")
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
# Manual processing from EDA, 
# know some of these parameters need to be changed
#...............................................................................................
txs$learnerlib[txs$target == "precip_mean_cont_scale_clst"] <- list(list(mlr::makeLearner("regr.lm", predict.type = "response")))
txs$learnerlib[txs$target == "wlthrcde_combscor_cont"] <- list(list(mlr::makeLearner("regr.lm", predict.type = "response")))
txs$learnerlib[txs$target == "ITN_fctb"] <- list(list(mlr::makeLearner("classif.logreg", predict.type = "prob")))


anemiamset<- c("hab57_fctb", "anyatm_cont_logit_scale_clst", "farmer_fctb",
                "hiv03_fctb", "hv009_cont_scale", "hv21345_fctb",
                "ITN_fctb", "wlthrcde_combscor_cont")

txs$task[txs$target == "hab57_fctb"] <- list(mlr::makeClassifTask(target = "hab57_fctb", data = dt.ml.cc[,anemiamset], positive = "no"))



farmermset <- c("farmer_fctb", "anyatm_cont_logit_scale_clst", "hab57_fctb",
                "hiv03_fctb", "hv009_cont_scale", "hv21345_fctb",
                "ITN_fctb", "wlthrcde_combscor_cont")

txs$task[txs$target == "farmer_fctb"] <- list(mlr::makeClassifTask(target = "farmer_fctb", 
                                                                   data = dt.ml.cc[,farmermset], positive = "farmer"))

eduset <- c("hv108_cont_scale", "urban_rura_fctb", "hv104_fctb")
txs$task[txs$target == "hv108_cont_scale"] <- list(mlr::makeRegrTask(target = "hv108_cont_scale", 
                                                                           data = dt.ml.cc[,eduset]))





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
                                                 time = "11-00:00:00"))

cat("*************************** \n Submitted SL models \n *************************** ")
























