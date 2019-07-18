#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make my own performance measure to estimate how well
# our IPTWeights are balancing the baseline covariates
# This is inspired by/is a slight extension of the function/methods presented in 
# Y. Zhu et al. "A Boosting Algorithm ..." 2015
#----------------------------------------------------------------------------------------------------
source("R/00-functions_epi.R")
library(mlr)
#..............................................................
# Binary Treatment Case
#..............................................................
my.covarbal.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame for Dij calculations
  dat <- mlr::getTaskData(task)
  n <- nrow(dat)
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  
  type <- mlr::getTaskType(task)
  if(type == "classif"){
    if( length(mlr::getTaskDesc(task)$class.levels) > 2){
      stop("This function does not support categorical treatment types")
    }
  }

  
  # find the type for the iptws from the task
 # warning("Type for IPTW Weight Calculations are being decided on class of target within data. Need to ensure proper coding")
  type <- ifelse(type == "classif", "binary",
                 ifelse(type == "regr", "continuous", NA))
  
  
  # pred <- mlr::getPredictionProbabilities(pred)
  
  # get inverse probability weights
  wi <- get_iptw_prob(task = task, preds = pred, type = type)
  
  # coerce back to vector
  nulldist <- unname( unlist(nulldist) )
  
  if(!is.atomic(nulldist)){
    stop("There is an issue with coercion from mlr performance to 
         a vector that we can call")
  }
  
  
  # pull apart data to find Djs
  data.list <- lapply(1:length(covars), 
                      function(x){return(as.data.frame(dat[, c(target, covars[x]) ]))})
  
  data.dist <- lapply(data.list, function(x){
    
    if(is.factor(x[,1])){
      x[,1] <- as.numeric(x[,1])
    }
    
    if(is.factor(x[,2])){
      x[,2] <- as.numeric(x[,2])
    }
    
    nrand <- floor(nrow(x) * 0.1) # take a 10% sample for speed
    nrand.rows <- sample(x = 1:nrow(x), size = nrand, replace = F) # pull random rows
    ret <- energy::dcor(x = x[nrand.rows, 1], y = x[nrand.rows, 2])
    
    return(ret)
  })
  
  # calculate average dist
  dij <- unlist(data.dist)
  dij <- mean(dij)
  
  # now calculate Z score
  Z <- ( dij - mean(nulldist) ) / sd(nulldist)
  
  return(Z)
  }



# Generate the Measure object for binary tx
my.covarbal = mlr::makeMeasure(
  id = "my.covarbal", 
  name = "Baseline covariate balance estimator",
  properties = c("classif", "regr"),
  extra.args = list(nulldist = NA),
  minimize = TRUE, 
  best = -Inf, worst = Inf,
  fun = my.covarbal.fun
)

