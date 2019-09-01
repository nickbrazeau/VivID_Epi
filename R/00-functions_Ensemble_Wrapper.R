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



# make the ensemble learner
#' @description 
#' simple approach with avergae as method of combining

make_avg_Stack <- function(
  task = task, 
  learners = learners.list){
  
  if(mlr::getTaskType(task) == "classif"){
    learners.list <- baselearners.list$classif
    baselearners <- lapply(learners.list, makeLearner, predict.type = "prob")
    m = makeStackedLearner(base.learners = baselearners,
                           predict.type = "prob",
                           method = "average")
    
    
  } else if(mlr::getTaskType(task) == "regr"){
    learners.list <- baselearners.list$regress
    
    baselearners <- lapply(learners.list, makeLearner, predict.type = "response")
    m <- makeStackedLearner(base.learners = baselearners,
                            predict.type = "response",
                            method = "average")
    
  } else {
    stop("You must have type binary or continuous for access to learners")
  }
  
  return(m)
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












  