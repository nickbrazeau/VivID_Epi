source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")
source("R/00-make_null_IPTW_distribs_brownian.R")
source("R/00-my_IPTW_performance_measure_energy.R")
set.seed(48, "L'Ecuyer")
library(tidyverse)
library(mlr)
library(rslurm)

dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
#........................
# manipulate tx map
#........................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(target = column_name,
                positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", target), "continuous",
                              ifelse(grepl("fctb", target), "binary", NA)))
  
#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
dt.ml <- dt %>% 
  dplyr::select(c("pv18s" , "pfldh", "hv005_wi", txs$target, 
                "longnum", "latnum"))

# subset to complete cases
dt.ml <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  data.frame(.)


#........................
# Subset and Store Dataframes for Tasks
#........................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml %>% 
                            dplyr::select(c(x, y, "longnum", "latnum"))
                          return(ret)
                        })

#........................
# Pull Out Coords
#........................
txs$coordinates <- purrr::map(txs$data, function(x){
  ret <- x %>% 
    dplyr::select(c("longnum", "latnum"))
  return(ret)
})

#-------------------------------------------------
# Setup Models
#-------------------------------------------------

# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type", "coordinates")], 
                         .f = make_class_task)

# # now look and correct class imbalance
# txs$task <- purrr::pmap(txs[,c("task", "type")], 
#                          .f = find_Class_Imbalance,
#                         classimb_tol = 0.6,
#                         smotenn = 5)


# now make the ensemble learner
txs$learner <- purrr::map(txs$task, make_simple_Stack, 
                          learners = baselearners.list)



###################################################
###################################################
######     Set HyperParameters for Tuning    ######     
###################################################
###################################################
# mlr::listLearners()
# ParamHelpers::getParamSet
# Note, issue with ksvm -- https://stackoverflow.com/questions/15895897/line-search-fails-in-training-ksvm-prob-model
#
# L1/L2 Regularization (glmnet), alpha: The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as (1-α)/2||β||_2^2+α||β||_1 alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
# K-Nearest Neighbors, k: Number of neighbors considered.
# Single Vector Machine, cost: cost of constraints violation (default: 1)—it is the ‘C’-constant of the regularization term in the Lagrange formulation.
# GAMBoost -- not going to tune 
# Random Forest, Mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 

# make a parameter set to explore
hyperparams_to_tune.regr <- ParamHelpers::makeParamSet(
  makeNumericParam("regr.glmnet.alpha", lower = 0, upper = 1),
  makeNumericParam("regr.kknn.k", lower = 1, upper = 26 ),
  makeNumericParam("regr.svm.cost", lower = 1, upper = 5),
  makeNumericParam("regr.randomForest.mtry", lower = 1, upper = 10 )
)

hyperparams_to_tune.regr.grid <-c(
  "regr.glmnet.alpha" = 11L,
  "regr.kknn.k" = 6L,
  "regr.svm.cost" = 5L,
  "regr.randomForest.mtry" =  10L 
 )


hyperparams_to_tune.classif <- ParamHelpers::makeParamSet(
  makeNumericParam("classif.glmnet.alpha", lower = 0, upper = 1),
  makeNumericParam("classif.kknn.k", lower = 1, upper = 26 ),
  makeNumericParam("classif.svm.cost", lower = 1, upper = 5),
  makeNumericParam("classif.randomForest.mtry", lower = 1, upper = 10 )
)



hyperparams_to_tune.classif.grid <-c(
  "classif.glmnet.alpha" = 11L,
  "classif.kknn.k" = 6L,
  "classif.svm.cost" = 5L,
  "classif.randomForest.mtry" =  10L 
)



txs$hyperparam <- purrr::map(txs$type, function(x){
  
  if(x == "continuous"){
    return(hyperparams_to_tune.regr)
  } else if(x == "binary"){
    return(hyperparams_to_tune.classif)
  } else {
    stop("Type not specified")
  }
  
})

# Make a Grid to Search On
txs$ctrl <- purrr::map(txs$type, function(x){
  
  if(x == "continuous"){
    return(makeTuneControlGrid(resolution = hyperparams_to_tune.regr.grid))
    
  } else if(x == "binary"){
    return(makeTuneControlGrid(resolution = hyperparams_to_tune.classif.grid))
    
  } else {
    stop("Type not specified")
  }
  
})


###################################################
###################################################
######   Set Performance Metrics/Null Dist   ######     
###################################################
###################################################
nulldist <- readRDS("analyses/06-IPW_ML/00-null_distributions/null_dist_return.RDS")

txs <- dplyr::left_join(txs, nulldist, by = "target")

txs$performmeasure <- lapply(1:nrow(txs), function(x) return(my.covarbal))

# set the null distribution for each respective DAG
txs$performmeasure <- map2(txs$performmeasure, txs$nulldist, function(x, y){
  ret <- mlr::setMeasurePars(x, 
                             par.vals = list(nulldist = y))
  })


###################################################
###################################################
######        RESAMPLING APPROACH            ######     
###################################################
###################################################
# resampling approach with spatial CV considered
# https://mlr.mlr-org.com/articles/tutorial/handling_of_spatial_data.html
# https://mlr.mlr-org.com/articles/tutorial/resample.html -- this is 3-fold CV (even though it is iters)
rdesc <- makeResampleDesc("SpCV", iters = 3)
txs$rdesc <- lapply(1:nrow(txs), function(x) return(rdesc))





#........................................................................
# SLURM 
#........................................................................

slurm_tunemodel <- function(learner, task, rdesc, hyperparam, ctrl, performmeasure){
  
  # setup parallelization
  target <- mlr::getTaskTargets(task)
  storagedir <- paste0(getwd(), "/", target)
  base::dir.create(path = storagedir, recursive = F) # this should only create one level
  
  parallelMap::parallelStartBatchtools(storagedir = storagedir, 
                                       bt.resources = list(walltime = 432000, ncpus = 8))
  
  ret <- mlr::tuneParams(learner = learner, 
                    task = task, 
                    resampling = rdesc, 
                    par.set = hyperparam,
                    control = ctrl,
                    measures = performmeasure, 
                    show.info = T)
  
  parallelMap::parallelStop()
  
  return(ret)
  
}

paramsdf <- txs[,c("learner", "task", "rdesc", "hyperparam", "ctrl", "performmeasure")]

# for slurm on LL
setwd("analyses/06-IPW_ML/tune_modelparams")
sjob <- rslurm::slurm_apply(f = slurm_tunemodel, 
                    params = paramsdf, 
                    jobname = 'vivid_preds',
                    nodes = 18, 
                    cpus_per_node = 8, 
                    submit = T,
                    slurm_options = list(mem = 128000,
                                         'cpus-per-task' = 8,
                                         error =  "%A_%a.err",
                                         output = "%A_%a.out",
                                         time = "11-00:00:00"))


cat("*************************** \n Submitted tuning models \n *************************** ")
  

