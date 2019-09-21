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

get_iptw_prob <- function(task, SLpreds){
  
  if(mlr::getTaskType(task) == "classif"){
    
    # pull details from mlr for numerator
    pos.class <- mlr::getTaskDesc(task)$positive
    target <- mlr::getTaskTargetNames(task) 

    ps <- SLpreds
    exposure <- mlr::getTaskData(task)[, target]
    pexp <- mean(exposure == pos.class)
    
    iptw_s <- ifelse(exposure == pos.class,
                     pexp/ps,
                     (1-pexp)/(1-ps))
    
    
  } else if(mlr::getTaskType(task) == "regr"){
    
    # following assumptions in Robbins 2000/Zhu 2015 PMC4749263
    # Code chunk follows Appendix A of Zhu 2015 PMC4749263
    # and Hernan Causal inference, Chapt 12 -- Program 12.4 
    # note, Hernan doesn't put on standard normal though
    preds <- SLpreds
    target <- mlr::getTaskTargetNames(task) 
    exposure <- mlr::getTaskData(task)[, target]
    model.num <- lm(exposure~1) 

    # divide residuals by model variance to put on standard normal
    ps.num <- dnorm(
      (exposure - model.num$fitted)/summary(model.num)$sigma,
      mean = 0, sd = 1, log = F)
    
    # get denominator 
    # but first am standardizing each value on to the 
    # standard normal for stabilization purposes
    ps.denom <- dnorm(
      (exposure - preds)/( sd( (exposure - preds)) ),
      mean = 0, sd = 1, log = F)
    
    iptw_s <- ps.num/ps.denom
    
  } else{
    stop("Type must be either binary or continous")
  }
  
  return(as.vector(iptw_s))
  
}


