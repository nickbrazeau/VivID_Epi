#----------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate find the best set of hyperparameters for our data
# Notes: Because I have done a grid.search approach, this optimization problem has 
# become "ridiculously" parallel. As a result, I have "overwritten" the default
# TuneControl objects to better parallelize the search on slurm 
# (notes to future self -- you tried to use `mlr`'s parallelization approach with parallelmap but could not get it to cooperate on slurm) 
#----------------------------------------------------------------------------------------------------


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

#--------------------------------------
# Setup tasks & stacked learner
#--------------------------------------

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



#--------------------------------------
# Setup tasks & stacked learner
#--------------------------------------
nulldist <- readRDS("analyses/06-IPW_ML/00-null_distributions/null_dist_return.RDS")

txs <- dplyr::left_join(txs, nulldist, by = "target")

txs$performmeasure <- lapply(1:nrow(txs), function(x) return(my.covarbal))

# set the null distribution for each respective DAG
txs$performmeasure <- map2(txs$performmeasure, txs$nulldist, function(x, y){
  ret <- mlr::setMeasurePars(x, 
                             par.vals = list(nulldist = y))
})


#--------------------------------------
# Setup resampling
#--------------------------------------
# resampling approach with spatial CV considered
# rdesc <- makeResampleDesc("SpRepCV", fold = 5, reps = 5)
rdesc <- makeResampleDesc("SpCV", iters = 3)
txs$rdesc <- lapply(1:nrow(txs), function(x) return(rdesc))

#--------------------------------------
# Setup dummy col
#--------------------------------------
txs$tune <- T



####################################################################################
####################################################################################
###########   Make a Hyperparamters for Tuning Submission Matrix     ###############     
####################################################################################
####################################################################################
# mlr::listLearners()
# ParamHelpers::getParamSet
#
# L1/L2 Regularization (glmnet), alpha: The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as (1-α)/2||β||_2^2+α||β||_1 alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
# K-Nearest Neighbors, k: Number of neighbors considered.
# Single Vector Machine, C: cost of constraints violation (default: 1) this is the `C'-constant of the regularization term in the Lagrange formulation.
# Single Vector Machine, kernel: The kernel function used in training and predicting. This parameter can be set to any function, of class kernel, which computes the inner product in feature space between two vector arguments (see kernels). kernlab provides the most popular kernel functions which can be used by setting the kernel parameter to the following strings: 
# GAMBoost -- not going to tune 
# TODO consider gamboost boosting for tuning
# Random Forest, Mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 

# make a parameter set to explore
hyperparams_to_tune <- ParamHelpers::makeParamSet(
  makeNumericParam("glmnet.alpha", lower = 0, upper = 1),
  makeNumericParam("kknn.k", lower = 1, upper = 30 ),
  makeNumericParam("svm.cost", lower = 1, upper = 5),
  makeNumericParam("randomForest.mtry", lower = 1, upper = 10 )
)

hyperparams_to_tune.ctrl <-c(
  "glmnet.alpha" = 2L, #11L
  "kknn.k" = 1L, #7L
  "svm.cost" = 1L, #5L
  "randomForest.mtry" =  1L #10L
)

hyperparams_to_tune.grid <- ParamHelpers::generateGridDesign(hyperparams_to_tune, 
                                                             hyperparams_to_tune.ctrl) %>% 
  dplyr::mutate(tune = T) %>% 
  dplyr::mutate(kknn.k = floor(kknn.k))



#--------------------------------------
# combine learners, tasks, resampling and hyperparams
#--------------------------------------

txs.hyperparams <- txs[, c("target", "learner", "task", "rdesc", "performmeasure", "tune")]

txs.hyperparams <- dplyr::left_join(txs.hyperparams, 
                                    hyperparams_to_tune.grid,
                                    by = "tune") %>% 
  dplyr::select(-c("tune"))



#--------------------------------------
# Rewrite Learners to Account for Hyperparams
#--------------------------------------

liftover_hyperparams <- function(learner, task,
                                 # now hyperparams
                                 glmnet.alpha,
                                 kknn.k,
                                 svm.cost,
                                 randomForest.mtry){
  
  # get task type to set learner
  type <- mlr::getTaskType(task)
  parvals <- list(glmnet.alpha, kknn.k, svm.cost, randomForest.mtry)
  names(parvals) <- paste0(type, c(".glmnet.alpha", ".kknn.k", ".svm.cost", ".randomForest.mtry"))
  learner.hp <- mlr::setHyperPars(learner, par.vals = parvals)
  return(learner.hp)
}

txs.hyperparams$learner <- purrr::pmap(txs.hyperparams[, c("learner", "task", 
                                                         "glmnet.alpha",
                                                         "kknn.k",
                                                         "svm.cost",
                                                         "randomForest.mtry")],
                                       liftover_hyperparams)

#--------------------------------------
# Now that we have modified learners,
# can subset to a smaller paramdf
#--------------------------------------
txs.hyperparams <- txs.hyperparams[, c("target", "learner", "task", "rdesc", "performmeasure")]


#........................................................................
# SLURM 
#........................................................................

slurm_tunemodel <- function(target, learner, task, rdesc, performmeasure){
  
  # Calculate the performance measures
  ret <- mlr::resample(learner = learner, 
                       task = task, 
                       resampling = rdesc, 
                       measures = performmeasure, 
                       show.info = T)
  
  return(ret)
  
}


# for slurm on LL
setwd("analyses/06-IPW_ML/tune_modelparams")
ntry <- 50
sjob <- rslurm::slurm_apply(f = slurm_tunemodel, 
                            params = txs.hyperparams, 
                            jobname = 'vivid_preds',
                            nodes = 1, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry - 1, 
                                                                 50),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))


cat("*************************** \n Submitted tuning models \n *************************** ")


