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

  
  
  # pred <- mlr::getPredictionProbabilities(pred)
  
  # get inverse probability weights 
  # wi <- get_iptw_prob(task = task, preds = pred, type = type) 
  # function having trouble on server bc of sub-nodes,
  # so had to copy and paste it here 
  get_iptw_prob <- function(task, preds){
    
    if(mlr::getTaskType(task) == "classif"){
      
      # pull details from mlr for numerator
      pos.class <- mlr::getTaskDesc(task)$positive
      target <- mlr::getTaskTargetNames(task) 
      data <- mlr::getTaskData(task)
      
      
      ps <- mlr::getPredictionProbabilities(preds)
      exposure <- mlr::getPredictionTruth(preds)
      pexp <- mean(exposure == pos.class)
      
      iptw_u <- ifelse(exposure == pos.class,
                       1/ps,
                       1/(1-ps))
      
      iptw_s <- pexp*iptw_u
      
    } else if(mlr::getTaskType(task) == "regr"){
      
      preds <- preds$data
      
      # following assumptions in Robbins 2000/Zhu 2015 PMC4749263
      model.num <- lm(preds$truth~1) # this is the intercept in the classic way
      #TODO check if sigma truly necessary 
      # divide residuals by model variance to put on standard normal
      ps.num <- dnorm(
        (preds$truth - model.num$fitted)/summary(model.num)$sigma,
        mean = 0, sd = 1, log = F)
      
      ps.denom <- dnorm(
        (preds$truth - preds$response)/( sd( (preds$truth - preds$response)) ),
        mean = 0, sd = 1, log = F)
      
      iptw_s <- ps.num/ps.denom
      
    } else{
      stop("Type must be either binary or continous")
    }
    
    return(iptw_s)
    
  } # copying and pasting function from R/00-functions_epi.R here
  wi <- get_iptw_prob(task = task, preds = pred)
  
  
  
  
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

