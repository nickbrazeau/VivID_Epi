# https://mlr.mlr-org.com/articles/tutorial/usecase_regression.html
library(tidyverse)
library(mlr)
data(BostonHousing2, package = "mlbench")
glimpse(BostonHousing2)


baselearners <- c("regr.lm",
                  "regr.glmnet", 
                  "regr.kknn",
                  "regr.ksvm",
                  "regr.randomForest")
covars <- c("medv", "crim", "zn", "indus", "chas", "nox", "rm", "age", "dis", "rad", "tax")

# task to fit
regr.task = makeRegrTask(data = BostonHousing2[,covars], target = "crim",
                         coordinates = BostonHousing2[,c("lon", "lat")])
regr.task

# make our stacked learner/ensemble learner
stck.lrnr <- makeStackedLearner(base.learners = baselearners,
                                predict.type = "response", 
                                method = "average")

mod = train(stck.lrnr, regr.task)
preds <- predict(mod, regr.task)

