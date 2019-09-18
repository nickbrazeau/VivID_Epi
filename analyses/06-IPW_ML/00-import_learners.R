baselearners.list <- list(
  classif =  c("classif.logreg",
               "classif.cvglmnet", 
               "classif.cvglmnet", 
               "classif.cvglmnet", 
               "classif.gamboost",
               "classif.svm",
               "classif.gausspr",
               "classif.ranger"),
  regress = c("regr.lm",
              "regr.cvglmnet", 
              "regr.cvglmnet", 
              "regr.cvglmnet", 
              "regr.svm",
              "regr.gamboost",
              "regr.gausspr",
              "regr.ranger")
)


base.learners.classif <- lapply(baselearners.list$classif, function(x) return(mlr::makeLearner(x, predict.type = "prob")))
base.learners.regr <- lapply(baselearners.list$regress, function(x) return(mlr::makeLearner(x, predict.type = "response")))


#.................
# manipulate glmnets classif
#.................
classifglmnets.num <- which(grepl("glmnet", baselearners.list$classif))
classifglmnets <-  base.learners.classif[ classifglmnets.num ]
classifglmnets[[1]] <- mlr::setHyperPars(classifglmnets[[1]], alpha = 1)
classifglmnets[[2]] <- mlr::setHyperPars(classifglmnets[[2]], alpha = 0.5)
classifglmnets[[3]] <- mlr::setHyperPars(classifglmnets[[3]], alpha = 0)

base.learners.classif[ classifglmnets.num ] <- classifglmnets

#.................
# manipulate glmnets classif
#.................
regrglmnets.num <- which(grepl("glmnet", baselearners.list$regress))
regrglmnets <-  base.learners.regr[ regrglmnets.num ]
regrglmnets[[1]] <- mlr::setHyperPars(regrglmnets[[1]], alpha = 1)
regrglmnets[[2]] <- mlr::setHyperPars(regrglmnets[[2]], alpha = 0.5)
regrglmnets[[3]] <- mlr::setHyperPars(regrglmnets[[3]], alpha = 0)

base.learners.regr[ regrglmnets.num ] <- regrglmnets



