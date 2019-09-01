#----------------------------------------------------------------------------------------------------
# Purpose of this script is to twofold:
# (1) Make iptw calculator for generalized propensity scores 
# (2) Make my own performance measure to estimate how well
#     our IPTWeights are balancing the baseline covariates
#    This is inspired by/is a slight extension of the function/methods presented in 
#     Y. Zhu et al. "A Boosting Algorithm ..." 2015
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R")
library(mlr)


#----------------------------------------------------------------------------------------------------
# iptw prob
#----------------------------------------------------------------------------------------------------

get_iptw_prob <- function(task, preds){
  
  if(mlr::getTaskType(task) == "classif"){
    
    # pull details from mlr for numerator
    pos.class <- mlr::getTaskDesc(task)$positive
    target <- mlr::getTaskTargetNames(task) 

    ps <- mlr::getPredictionProbabilities(preds, cl = pos.class)
    exposure <- mlr::getPredictionTruth(preds)
    pexp <- mean(exposure == pos.class)
    
    iptw_s <- ifelse(exposure == pos.class,
                     pexp/ps,
                     (1-pexp)/(1-ps))
    
    
  } else if(mlr::getTaskType(task) == "regr"){
    
    # following assumptions in Robbins 2000/Zhu 2015 PMC4749263
    # Code chunk follows Appendix A of Zhu 2015 PMC4749263
    # and Hernan Causal inference, Chapt 12 -- Program 12.4 
    # note, Hernan doesn't put on standard normal though
    preds <- preds$data
    model.num <- lm(preds$truth~1) 

    # divide residuals by model variance to put on standard normal
    ps.num <- dnorm(
      (preds$truth - model.num$fitted)/summary(model.num)$sigma,
      mean = 0, sd = 1, log = F)
    
    # get denominator 
    # but first am standardizing each value on to the 
    # standard normal for stabilization purposes
    ps.denom <- dnorm(
      (preds$truth - preds$response)/( sd( (preds$truth - preds$response)) ),
      mean = 0, sd = 1, log = F)
    
    iptw_s <- ps.num/ps.denom
    
  } else{
    stop("Type must be either binary or continous")
  }
  
  return(iptw_s)
  
}




#..............................................................
# My Covariance Function to assess how well our 
# calculated weights actually result in no baseline confounding
#..............................................................
my.covarbal.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame for Dij calculations
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  
  type <- mlr::getTaskType(task)
  if(type == "classif"){
    if( length(mlr::getTaskDesc(task)$class.levels) > 2){
      stop("This function does not support categorical treatment types")
    }
  }

  
  #........................
  # Get and Apply IPTWs
  #........................
  # Note, predications can have fewer observations than task
  # as this is part of the folds/iterations

  wi <- get_iptw_prob(task = task, preds = pred)
  
  dat <- mlr::getTaskData(task)[pred$data$id, ] # need to pull only those observations we predicted, so we can apply weights
  dat.weighted.rows <- sample(x = 1:nrow(dat), size = nrow(dat), prob = wi, replace = T)
  dat.weighted <- dat[dat.weighted.rows, ]
  

  #........................
  # Nulldist to useful by
  # coercing back to vector
  #........................
  
  nulldist <- unname( unlist(nulldist) )
  
  if(!is.atomic(nulldist)){
    stop("There is an issue with coercion from mlr performance to 
         a vector that we can call")
  }
  
  
  # pull apart data to find Djs
  data.list <- lapply(1:length(covars), 
                      function(x){return(as.data.frame(dat.weighted[, c(target, covars[x]) ]))})
  
  data.dist <- lapply(data.list, function(x){
    
    if(is.factor(x[,1])){ # note, only have binary factors so this ok
      x[,1] <- as.numeric(x[,1])
    }
    
    if(is.factor(x[,2])){
      x[,2] <- as.numeric(x[,2])
    }
    
    ret <- energy::dcor(x = x[, 1], y = x[, 2])
    
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

