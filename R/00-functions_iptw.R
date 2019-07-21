#----------------------------------------------------------------------------------------------------
# iptw prob
#----------------------------------------------------------------------------------------------------

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
  
}



#----------------------------------------------------------------------------------------------------
# tune our learner 
#----------------------------------------------------------------------------------------------------

tune_stacked_learner <- function(learner, task, tuneresult){
  if(mlr::getTaskType(task) == "classif"){
    stck.lrnr.tuned <- setHyperPars(learner, 
                                    classif.glmnet.alpha = tuneresult$x$classif.glmnet.alpha,
                                    classif.kknn.k = tuneresult$x$classif.kknn.k,
                                    classif.ksvm.C = tuneresult$x$classif.ksvm.C,
                                    classif.ksvm.kernel = tuneresult$x$classif.ksvm.kernel,
                                    classif.randomForest.mtry = tuneresult$x$classif.randomForest.mtry
    )
  } else if(mlr::getTaskType(task) == "regr"){
    
    stck.lrnr.tuned <- setHyperPars(learner, 
                                    regr.glmnet.alpha = tuneresult$x$regr.glmnet.alpha,
                                    regr.kknn.k = tuneresult$x$regr.kknn.k,
                                    regr.ksvm.C = tuneresult$x$regr.ksvm.C,
                                    regr.ksvm.kernel = tuneresult$x$regr.ksvm.kernel,
                                    regr.randomForest.mtry = tuneresult$x$regr.randomForest.mtry
    )
  }
  stck.lrnr.tuned <- setLearnerId(stck.lrnr.tuned, "stacked_learner_tuned")
  return(stck.lrnr.tuned)
}









