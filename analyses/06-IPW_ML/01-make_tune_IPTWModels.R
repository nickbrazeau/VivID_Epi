source("R/00-simple_Ensemble_Wrapper.R")
source("R/00-functions_epi.R")
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
txs$learner <- purrr::map(txs$type, make_simple_Stack, 
                          learners = baselearners.list)



###################################################
###################################################
######     Set HyperParameters for Tuning    ######     
###################################################
###################################################

# L1/L2 Regularization (glmnet), alpha: The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as (1-α)/2||β||_2^2+α||β||_1 alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
# K-Nearest Neighbors, k: Number of neighbors considered.
# Single Vector Machine, C: cost of constraints violation (default: 1) this is the `C'-constant of the regularization term in the Lagrange formulation.
# Single Vector Machine, kernel: The kernel function used in training and predicting. This parameter can be set to any function, of class kernel, which computes the inner product in feature space between two vector arguments (see kernels). kernlab provides the most popular kernel functions which can be used by setting the kernel parameter to the following strings: 
# GAMBoost -- not going to tune 
# TODO consider gamboost boosting for tuning
# Random Forest, Mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 

# make a parameter set to explore
hyperparams_to_tune <- ParamHelpers::makeParamSet(
  makeNumericParam("regr.glmnet.alpha", lower = 0, upper = 1),
  makeNumericParam("regr.kknn.k", lower = 1, upper = 15 ),
  makeNumericParam("regr.ksvm.C", lower = 1, upper = 10),
  makeDiscreteParam("regr.ksvm.kernel", values = c("rbfdot", "polydot", "vanilladot", "besseldot", "splinedot", "stringdot")),
  makeNumericParam("regr.randomForest.mtry", lower = 1, upper = 10 )
)

txs$hyperparam <- lapply(1:nrow(txs), function(x) return(hyperparams_to_tune))

# Make a Grid to Search On
ctrl <- makeTuneControlGrid(resolution = 10)
txs$ctrl <- lapply(1:nrow(txs), function(x) return(ctrl))



###################################################
###################################################
######   Set Performance Metrics/Null Dist   ######     
###################################################
###################################################
nulldist <- readRDS("analyses/06-IPW_ML/00-null_distributions/null_dist_return.RDS")

txs <- dplyr::left_join(txs, nulldist, by = "target")

txs$performmeasure <- lapply(1:nrow(txs), function(x) return(my.covarbal))

txs$performmeasure <- map2(txs$performmeasure, txs$nulldist, function(x, y){
  # set the null distribution for each respective DAG
  ret <- mlr::setMeasurePars(x, 
                             par.vals = list(nulldist = y))
  })


###################################################
###################################################
######        RESAMPLING APPROACH            ######     
###################################################
###################################################
# resampling approach with spatial CV considered
rdesc <-makeResampleDesc("SpRepCV", fold = 5, reps = 5)
txs$rdesc <- lapply(1:nrow(txs), function(x) return(rdesc))





#........................................................................
# SLURM 
#........................................................................

slurm_tunemodel <- function(learner, task, rdesc, hyperparam, ctrl, performmeasure){
  
  ret <- mlr::tuneParams(learner = learner, 
                    task = task, 
                    resampling = rdesc, 
                    par.set = hyperparam,
                    control = ctrl,
                    measures = performmeasure, 
                    show.info = T)
  
  return(ret)
  
}

paramsdf <- txs[,c("learner", "task", "rdesc", "hyperparam", "ctrl", "performmeasure")]

# for slurm on LL
setwd("analyses/06-IPW_ML/")
ntry <- 18
sjob <- rslurm::slurm_apply(f = slurm_tunemodel, 
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
                                         time = "5-00:00:00"))

# saveRDS(object = txs, file = "results/ensemble_predictive_probs.RDS")





