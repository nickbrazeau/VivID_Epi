baselearners.list <- list(
  classif =  c("classif.logreg",
  #             "classif.glmnet", 
               "classif.gamboost",
               "classif.kknn",
               "classif.svm",
  #             "classif.gausspr",
               "classif.randomForest"),
  regress = c("regr.lm",
   #           "regr.glmnet", 
              "regr.gamboost",
              "regr.kknn",
              "regr.svm",
   #           "regr.gausspr",
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

#' @param numeric; The rate of oversampling is set internally based on prop imbalance
#' via 1 + max(props) - min(props)
#' @description  find class imbalance and adjust the task
find_Class_Imbalance <- function(task, type, classimb_tol = 0.6,
                                 smotenn = 5){
  if(type == "binary"){
    
    props <- prop.table(table(mlr::getTaskTargets(task)))
    
    if(max(props) > classimb_tol){
      message("There appears to be class-imbalance for ", task$task.desc$target)
      smoterate = 1 + max(props) - min(props)
      
      task <- mlr::smote(task, rate = smoterate, nn = smotenn)
      return(task)
      
    } else {
      return(task)
    } # end inner ifelse
    
  } else {
    return(task)
  } # end outer ifelse
}





# make the ensemble learner
#' @description 
#' simple approach with hill climb as method of combining

make_hillclimb_Stack <- function(
  task = task, 
  learners = learners.list){
  
  if(mlr::getTaskType(task) == "classif"){
    learners.list <- baselearners.list$classif
    baselearners <- lapply(learners.list, makeLearner, predict.type = "prob")
    m = makeStackedLearner(base.learners = baselearners,
                           predict.type = "prob",
                           method = "hill.climb")
    
    
  } else if(mlr::getTaskType(task) == "regr"){
    learners.list <- baselearners.list$regress
    
    baselearners <- lapply(learners.list, makeLearner, predict.type = "response")
    m <- makeStackedLearner(base.learners = baselearners,
                            predict.type = "response",
                            method = "hill.climb")
    
  } else {
    stop("You must have type binary or continuous for access to learners")
  }
  
  return(m)
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




  