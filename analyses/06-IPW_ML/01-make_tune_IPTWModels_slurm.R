#----------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate find the best set of hyperparameters for our data
# Notes: Here I am using mlr built-in multiplex option and individual tuning my learners 
# will take the intersection of the best hyperparameter set and apply it to our stacked ensemble
#----------------------------------------------------------------------------------------------------

source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-IPTW_functions.R")
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
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA))) %>% 
  dplyr::rename(target = column_name)


#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
varstoinclude <- c("pv18s" , "pfldh", "hv005_wi", txs$target,
                   "alt_dem_cont_scale_clst", "hab1_cont_scale", "hv104_fctb", "wtrdist_cont_log_scale_clst", # need to add in covariates that don't have confounding ancestors but are needed elsewhere
                   "hiv03_fctb", # no longer considered risk factor bc too few observations
                   "longnum", "latnum")
dt.ml <- dt %>% 
  dplyr::select(varstoinclude)

# subset to complete cases
dt.ml.cc <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum")) %>% 
  data.frame(.)

dt.ml.coords <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  data.frame(.)

#........................
# Subset Data for Tuning
# 80% of data to train on 
# to avoid overfitting
#........................
nrows.df <- nrow(dt.ml.cc)
pull <- sort( sample(1:nrows.df, size = 0.8*nrows.df) )

dt.ml.cc <- dt.ml.cc[pull, ]
dt.ml.coords <- dt.ml.coords[pull, ]

#........................
# Subset and Store Dataframes for Tasks
#........................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml.cc %>% 
                            dplyr::select(c(x, y))
                          return(ret)
                        })

#........................
# Pull Out Coords
#........................
txs$coordinates <- lapply(1:nrow(txs), function(x) return(dt.ml.coords))


#--------------------------------------
# Setup tasks & base learners
#--------------------------------------
# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type", "coordinates")], 
                        .f = make_class_task)



# now make the individual multiplexed base learners
base.learners.regr <- lapply(baselearners.list$regress, function(x) return(mlr::makeLearner(x, predict.type = "response")))
base.learners.regr <- mlr::makeModelMultiplexer(base.learners.regr)

base.learners.classif <- lapply(baselearners.list$classif, function(x) return(mlr::makeLearner(x, predict.type = "prob")))
base.learners.classif <- mlr::makeModelMultiplexer(base.learners.classif)

txs$learner <- purrr::map(txs$type, function(x){
    if(x == "continuous"){
      return(base.learners.regr)
    } else if (x == "binary"){
      return(base.learners.classif)
    }
  })



#--------------------------------------
# Setup null distributions
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
rdesc <- makeResampleDesc("SpRepCV", fold = 5, reps = 5)
txs$rdesc <- lapply(1:nrow(txs), function(x) return(rdesc))



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
# Random Forest, Mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 

# make a parameter set to explore
# REGRESSSION
hyperparams_to_tune.regr <- mlr::makeModelMultiplexerParamSet(
  base.learners.regr, 
  makeNumericParam("alpha", lower = 0, upper = 1),
  makeNumericParam("k", lower = 2, upper = 30 ), # knn of 1 just memorizes data basically
  makeNumericParam("cost", lower = 1, upper = 5),
  makeNumericParam("mtry", lower = 1, upper = 5 )
)

hyperparams_to_tune.regr.ctrl <-c(
  "regr.glmnet.alpha" = 11L, #11L
  "regr.kknn.k" = 29L, #7L
  "regr.svm.cost" = 5L, #5L
  "regr.randomForest.mtry" =  5L #10L
)

hyperparams_to_tune.regr.ctrl <- mlr::makeTuneControlDesign(design = 
  ParamHelpers::generateGridDesign(par.set = hyperparams_to_tune.regr,
                                   resolution = hyperparams_to_tune.regr.ctrl)
  )

# CLASSIFICATION
hyperparams_to_tune.classif <- mlr::makeModelMultiplexerParamSet(
  base.learners.classif, 
  makeNumericParam("alpha", lower = 0, upper = 1),
  makeNumericParam("k", lower = 2, upper = 30 ), # knn of 1 just memorizes data basically
  makeNumericParam("cost", lower = 1, upper = 5),
  makeNumericParam("mtry", lower = 1, upper = 5 )
)

hyperparams_to_tune.classif.ctrl <-c(
  "classif.glmnet.alpha" = 11L, #11L
  "classif.kknn.k" = 29L, #29L
  "classif.svm.cost" = 5L, #5L
  "classif.randomForest.mtry" =  5L #10L
)


hyperparams_to_tune.classif.ctrl <- mlr::makeTuneControlDesign(design = 
  ParamHelpers::generateGridDesign(par.set = hyperparams_to_tune.classif,
                                   resolution = hyperparams_to_tune.classif.ctrl)
  )

#--------------------------------------
#add in hyperparams and control   
#--------------------------------------
txs$hyperparams <- purrr::map(txs$type, function(x){
  if(x == "continuous"){
    return(hyperparams_to_tune.regr)
  } else if(x == "binary"){
    return(hyperparams_to_tune.classif)
  }
})

txs$ctrl <- purrr::map(txs$type, function(x){
  if(x == "continuous"){
    return(hyperparams_to_tune.regr.ctrl)
  } else if(x == "binary"){
    return(hyperparams_to_tune.classif.ctrl)
  }
})




#........................................................................
# SLURM 
#........................................................................

slurm_tunemodel <- function(learner, task, rdesc, hyperparams, ctrl, performmeasure){
  
  cat(mlr::getTaskTargetNames(task))
  ret <- mlr::tuneParams(learner = learner, 
                         task = task, 
                         resampling = rdesc, 
                         par.set = hyperparams,
                         control = ctrl,
                         measures = performmeasure, 
                         show.info = T)
  
  return(ret)
  
}

paramsdf <- txs %>% 
  dplyr::select(c("learner", "task", "rdesc", "hyperparams", "ctrl", "performmeasure")) 


# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = slurm_tunemodel, 
                            params = paramsdf, 
                            jobname = 'vivid_tunes_train',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry-1, 
                                                                 ntry),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))

cat(" *************************** \n Submitted tuning models \n *************************** ")


