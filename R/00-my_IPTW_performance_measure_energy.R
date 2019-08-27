#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make my own performance measure to estimate how well
# our IPTWeights are balancing the baseline covariates
# This is inspired by/is a slight extension of the function/methods presented in 
# Y. Zhu et al. "A Boosting Algorithm ..." 2015
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R")
library(mlr)
#..............................................................
# Binary Treatment Case
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

  
  
  # pred <- mlr::getPredictionProbabilities(pred)
  
  # get inverse probability weights 
  # adjusted as predictions have fewer observations than task
  # as this is part of the folds/iterations
  get_iptw_prob <- function(task, preds){
    
    if(mlr::getTaskType(task) == "classif"){
      
      # pull details from mlr for numerator
      pos.class <- mlr::getTaskDesc(task)$positive
      target <- mlr::getTaskTargetNames(task) 

      
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
  
  
  #........................
  # Apply IPTWs
  #........................
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
  
  #........................
  # CONSIDER DISTRIBUTION of IPTWs
  #........................
  
  
  # weights are on scale of 0, Inf, can log tranform these 
  # IMO, ideal distribution of weights would be ~ N(1, 0.25)
  # assumption of independence is violated here but we aren't using this for inference, just guidance
  # iptw.sc <- my.scale(iptw)
  # smsmpl <- sample(1:length(iptw.sc), 100)
  # ll <- sum(dnorm(x = iptw.sc[smsmpl], mean = 0, sd = 1, log = T))
  # too much data overwhelming likelihood... 
  # bc we are assuming 15000 iid draws...that are not iid
  
  iptw.fit <- MASS::fitdistr(wi, "normal")
  pd <- iptw.fit$estimate["mean"]/iptw.fit$estimate["sd"]
  
  # pull it in right direction
  if(Z > 1){
    z.pd <- Z/pd
  } else if(Z < 1) {
    z.pd <- Z * pd
  }
  
  
  return(z.pd)
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

