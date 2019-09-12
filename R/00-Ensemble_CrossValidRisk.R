


#' @details follows https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4262745/
#' 

ensemble_crossval_risk_pred <- function(learnerlib, task, proptrainset){
  
  #'
  #' @details Function to extract predictions as a list
  #' from the various individual learners
  
  get_preds <- function(trained, task, subset){
    if(mlr::getTaskType(task) == "classif"){
      
      # pull details from mlr 
      pos.class <- mlr::getTaskDesc(task)$positive
      preds.i <- predict(trained, task = task, subset = subset)
      preds.i <- mlr::getPredictionProbabilities(pred = preds.i, cl = pos.class)
      
    } else if(mlr::getTaskType(task) == "regr"){
      preds.i <- predict(trained, task = task, subset = subset)$data$response
    }
    return(preds.i)
  }
  
  
  #'
  #' @details This is the Cross validated risk -- function is essentially 
  #' lifted directly from SuperLearner::method.NNLS. I dropped observation weights
  #' and cahnge the return a bit
  nnls_cvrisk <- function (Z, Y, algnames, verbose = T){
    cvRisk <- apply(Z, 2, function(x){return( mean((x - Y)^2) )})
    
    names(cvRisk) <- algnames
    fit.nnls <- nnls::nnls( Z, Y)
    if (verbose) {
      message(paste("Non-Negative least squares convergence:", 
                    fit.nnls$mode == 1))
    }
    initCoef <- coef(fit.nnls)
    initCoef[is.na(initCoef)] <- 0
    if (sum(initCoef) > 0) {
      coef <- initCoef/sum(initCoef)
    }
    else {
      warning("All algorithms have zero weight", call. = FALSE)
      coef <- initCoef
    }
    out <- list(cvrisk = cvRisk, coef = coef)
    return(out)
  }
  
  
  
  #..............................
  # Setup pieces
  #..............................
  fulldat <- mlr::getTaskData(task)
  n <- mlr::getTaskSize(task)
  trainset <- sample(n, size = floor(proptrainset*n))
  
  # train base libraries
  learnerlib.trained <- purrr::map(learnerlib,
                                   function(x, task, trainset){
                                     ret <- mlr::train(learner = x, task = task, subset = trainset)
                                   }, task = task, trainset = trainset)
  # get validation set
  fulln <- seq(1, n, 1)
  valset <- fulln[! fulln %in% trainset]
  
  
  #..............................
  # prediction matrix, Zed
  #..............................
  Z <- matrix(NA, nrow = length(valset), ncol = length(learnerlib.trained))
  for(j in 1:length(learnerlib.trained)){
    Z[, j] <- get_preds(trained = learnerlib.trained[[j]],
                        task = task,
                        subset = valset)
  }
  
  #..............................
  # minimize cross validated risk
  #..............................
  algnames <- unlist( purrr::map(learnerlib, "name") )
  Y <- fulldat[valset, mlr::getTaskTargetNames(task)] # or really A
  if(is.factor(Y)){
    pos.class <- mlr::getTaskDesc(task)$positive
    Y <- ifelse(Y == pos.class, 1, 0)
  }
  
  if(length(learnerlib) == 1){ # can't do nnls on one alg
    cvrisk <- list(cvrisk = 1, coef = 1)
  } else{
    cvrisk <- nnls_cvrisk(Z = Z, Y = Y, algnames = algnames, verbose = T)
  }
  
  
  #..............................
  # Predict on Full Dataset
  #..............................
  # not full learner list makes it to final set
  
  finalalgspass <- cvrisk$coef > 0
  if(all(cvrisk$coef == 0)){stop("The CV Risk could not be calculated for this model")}
  
  
  learnerlib.trained.final <- learnerlib.trained[which(finalalgspass)]
  
  Zprime <- matrix(NA, nrow = nrow(fulldat), ncol = length(learnerlib.trained.final))
  for(j in 1:length(learnerlib.trained.final)){
    Zprime[, j] <- get_preds(trained = learnerlib.trained.final[[j]],
                             task = task,
                             subset = NULL
    )
  }
  
  EL.cvrisk.preds <- Zprime %*% cvrisk$coef[finalalgspass] 
  
  
  
  #..............................
  # outs
  #..............................
  ret <- list(cvrisk.coef = cvrisk$coef,
              alg.cvrisk.validationset = cvrisk$cvrisk,
              EL.predictions = unlist(EL.cvrisk.preds),
              task = task,
              Z = Z,
              Zprime = Zprime)
  
  return(ret)
  
  
}




