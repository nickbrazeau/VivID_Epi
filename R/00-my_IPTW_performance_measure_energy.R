#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make my own performance measure to estimate how well
# our IPTWeights are balancing the baseline covariates
# This is inspired by/is a slight extension of the function/methods presented in 
# Y. Zhu et al. "A Boosting Algorithm ..." 2015
#----------------------------------------------------------------------------------------------------
source("R/00-functions_epi.R")

#..............................................................
# Binary Treatment Case
#..............................................................
my.covarbal.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame for Dij calculations
  dat <- mlr::getTaskData(task)
  n <- nrow(dat)
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  # get inverse probability weights
  wi <- get_iptw_prob(task = task, preds = pred, type = "binary")
  
  # coerce back to vector
  nulldist <- unname( unlist(nulldist) )
  
  if(!is.atomic(nulldist)){
    stop("There is an issue with coercion from mlr performance to 
         a vector that we can call")
  }
  
  
  # pull apart data to find Djs
  data.list <- lapply(1:length(covars), 
                      function(x){return(as.data.frame(data[, c(target, covars[x]) ]))})
  
  data.dist <- lapply(data.list, function(x){
    
    if(is.factor(x[,1])){
      x[,1] <- as.numeric(x[,1])
    }
    
    if(is.factor(x[,2])){
      x[,2] <- as.numeric(x[,2])
    }
    
    ret <- energy::dcor(x = x[,1], y = x[,2])
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
my.covarbal = makeMeasure(
  id = "my.covarbal", 
  name = "Baseline covariate balance estimator",
  properties = c("classif", "regr"),
  extra.args = list(nulldist = NA),
  minimize = TRUE, 
  best = -Inf, worst = Inf,
  fun = my.covarbal.fun
)

