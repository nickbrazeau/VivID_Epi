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
    # same as mean(preds$truth)
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
findbesttuneresult <- function(path){
  res <- readRDS(path)
  dfres <- as.data.frame(res[[1]]$opt.path)
  dfres <- dfres %>% 
    dplyr::select(-c("dob", "eol", "error.message", "exec.time", "selected.learner")) %>% 
    tidyr::gather(., key = "hyperpar", val = "hyperparval", 1:(ncol(.)-1)) %>% 
    dplyr::filter(!is.na(hyperparval)) %>% 
    dplyr::group_by(hyperpar) %>% 
    dplyr::filter(my.covarbal.test.mean == min(my.covarbal.test.mean))
  
  return(dfres)
}



findbesttuneresult.simple <- function(path){
  res <- readRDS(path)
  dfres <- as.data.frame(res[[1]]$opt.path)
  
  simpletestmean <- colnames(dfres)[grepl("test.mean", colnames(dfres))]
  
  if(simpletestmean == "auc.test.mean"){
    dfres <- dfres %>% 
      dplyr::select(-c("dob", "eol", "error.message", "exec.time", "selected.learner")) %>% 
      tidyr::gather(., key = "hyperpar", val = "hyperparval", 1:(ncol(.)-1)) %>% 
      dplyr::filter(!is.na(hyperparval)) %>% 
      dplyr::group_by(hyperpar) %>% 
      dplyr::filter(auc.test.mean == min(auc.test.mean))
  } else if(simpletestmean == "mse.test.mean"){
    dfres <- dfres %>% 
      dplyr::select(-c("dob", "eol", "error.message", "exec.time", "selected.learner")) %>% 
      tidyr::gather(., key = "hyperpar", val = "hyperparval", 1:(ncol(.)-1)) %>% 
      dplyr::filter(!is.na(hyperparval)) %>% 
      dplyr::group_by(hyperpar) %>% 
      dplyr::filter(mse.test.mean == min(mse.test.mean))
  }
  
  # error catch
  # fit lowest model which is the first
  dfres <- dfres[!duplicated(dfres$hyperpar), ]
  
  return(dfres)
}




tune_stacked_learner <- function(learner, task, tuneresult){
  if(mlr::getTaskType(task) == "classif"){
    stck.lrnr.tuned <- setHyperPars(learner, 
                                    classif.glmnet.alpha = tuneresult$hyperparval[tuneresult$hyperpar == "classif.glmnet.alpha"],
                                    classif.kknn.k = tuneresult$hyperparval[tuneresult$hyperpar == "classif.kknn.k"],
                                    classif.svm.cost = tuneresult$hyperparval[tuneresult$hyperpar == "classif.svm.cost"],
                                    classif.randomForest.mtry = tuneresult$hyperparval[tuneresult$hyperpar == "classif.randomForest.mtry"]
    )
  } else if(mlr::getTaskType(task) == "regr"){
    
    stck.lrnr.tuned <- setHyperPars(learner, 
                                    regr.glmnet.alpha = tuneresult$hyperparval[tuneresult$hyperpar == "regr.glmnet.alpha"],
                                    regr.kknn.k = tuneresult$hyperparval[tuneresult$hyperpar == "regr.kknn.k"],
                                    regr.svm.cost = tuneresult$hyperparval[tuneresult$hyperpar == "regr.svm.cost"],
                                    regr.randomForest.mtry = tuneresult$hyperparval[tuneresult$hyperpar == "regr.randomForest.mtry"]
    )
  }
  stck.lrnr.tuned <- setLearnerId(stck.lrnr.tuned, "stacked_learner_tuned")
  return(stck.lrnr.tuned)
}









