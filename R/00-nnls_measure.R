my.nnls.fun <- function(task, model, pred, feats, extra.args){
  
  truth <- mlr::getPredictionTruth(pred)
  
  if(mlr::getTaskDesc(task)$type == "classif"){
    pos.class <- mlr::getTaskDesc(task)$positive
    truth <- ifelse(truth == pos.class, 1, 0)
    
  }
  

  
  if(mlr::getTaskDesc(task)$type == "regr") {
    Y <- mlr::getPredictionResponse(pred)
  } else if(mlr::getTaskDesc(task)$type == "classif") {
    Y <- mlr::getPredictionProbabilities(pred)
  }
  
  ret <- sum((truth - Y)^2)
  
  return(ret)

}



# Generate the Measure object for binary tx
my.nnls = mlr::makeMeasure(
  id = "my.covarbal", 
  name = "Baseline covariate balance estimator",
  properties = c("classif", "regr"),
  minimize = TRUE, 
  best = 0, worst = Inf,
  fun = my.nnls.fun
)
