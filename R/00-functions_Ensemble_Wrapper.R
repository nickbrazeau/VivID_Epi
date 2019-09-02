baselearners.list <- list(
  classif =  c("classif.logreg",
               "classif.glmnet", 
               "classif.gamboost",
               "classif.kknn",
               "classif.svm",
               "classif.gausspr",
               "classif.randomForest"),
  regress = c("regr.lm",
              "regr.glmnet", 
              "regr.gamboost",
              "regr.kknn",
              "regr.svm",
              "regr.gausspr",
              "regr.randomForest")
)




# make a task
make_class_task <- function(data = data, 
                            type = type,
                            target = target,
                            positive = positive,
                            coordinates = NULL){
  
  if(type == "binary"){
    task <- mlr::makeClassifTask(data = data, 
                                 target = target,
                                 positive = positive,
                                 coordinates = coordinates)
  } else if(type == "continuous"){
    task <- mlr::makeRegrTask(data = data, 
                              target = target,
                              coordinates = coordinates)
  }
  
  return(task)
  
}




#----------------------------------------------------------------------------------------------------
# Finding Tuning Result
#----------------------------------------------------------------------------------------------------
findbesttuneresult <- function(path){
  res <- readRDS(path)
  dfres <- as.data.frame(res[[1]]$opt.path)
  
  simpletestmean <- colnames(dfres)[grepl("test.mean", colnames(dfres))]
  
  if(simpletestmean == "logloss.test.mean"){
    dfres <- dfres %>% 
      dplyr::select(-c("dob", "eol", "error.message", "exec.time", "selected.learner")) %>% 
      tidyr::gather(., key = "hyperpar", val = "hyperparval", 1:(ncol(.)-1)) %>% 
      dplyr::filter(!is.na(hyperparval)) %>% 
      dplyr::group_by(hyperpar) %>% 
      dplyr::filter(logloss.test.mean == min(logloss.test.mean))
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


#----------------------------------------------------------------------------------------------------
# Apply Tuning Result
#----------------------------------------------------------------------------------------------------


tune_learner_library <- function(learnerlib, hyperparamstuned){
  
  lrnrnames <- learnerlib %>% 
    purrr::map(., "id") %>% 
    unlist(.)
  
  if(all(grepl("regr", lrnrnames))){
    
    learnerlib[[which(grepl("glmnet", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("glmnet", lrnrnames))]],
                                                                    alpha = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "regr.glmnet.alpha"],
                                                                    s = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "regr.glmnet.s"])
    
    learnerlib[[which(grepl("knn", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("knn", lrnrnames))]],
                                                                 k = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "regr.knn.k"])
    
    learnerlib[[which(grepl("svm", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("svm", lrnrnames))]],
                                                                 cost = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "regr.svm.cost"])
    
    learnerlib[[which(grepl("randomForest", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("randomForest", lrnrnames))]],
                                                                          mtry = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "regr.randomForest.mtry"])
    
  } else if(all(grepl("classif", lrnrnames))){
    learnerlib[[which(grepl("glmnet", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("glmnet", lrnrnames))]],
                                                                    alpha = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "classif.glmnet.alpha"],
                                                                    s = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "classif.glmnet.s"])
    
    learnerlib[[which(grepl("knn", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("knn", lrnrnames))]],
                                                                 k = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "classif.knn.k"])
    
    learnerlib[[which(grepl("svm", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("svm", lrnrnames))]],
                                                                 cost = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "classif.svm.cost"])
    
    learnerlib[[which(grepl("randomForest", lrnrnames))]] <- setHyperPars(learnerlib[[which(grepl("randomForest", lrnrnames))]],
                                                                          mtry = hyperparamstuned$hyperparval[hyperparamstuned$hyperpar == "classif.randomForest.mtry"])
  } else{
    stop("Not a task by mlr?")
  }
  
  return(learnerlib)
}















  