# https://mlr.mlr-org.com/articles/tutorial/usecase_regression.html

library(mlr)
data(BostonHousing2, package = "mlbench")
summary(BostonHousing2)
BostonHousing

baselearners <- c("regr.lm",
              "regr.glmnet",
              "regr.randomForest")

BostonHousing2 %>% 
  dplyr::group_by(town) %>%  # indus -- proportion of residential land for industry by town doesn't change -- like weather would not --> this isn't data leakage
  dplyr::summarise(
    meannox = mean(indus),
    sdnox = sd(indus)
  )

# ^ good sanity check for our data

colnames(BostonHousing2)
covar <- c("medv", "crim", "zn", "indus", "chas", "nox", "rm", "age", "dis", "rad", "tax")


# task to fit
regr.task = makeRegrTask(data = BostonHousing2[,covar], target = "crim",
                         coordinates = BostonHousing2[,c("lon", "lat")])
regr.task

# make our stacked learner/ensemble learner
stck.lrnr <- makeStackedLearner(base.learners = baselearners,
                                predict.type = "response", 
                                method = "average")
stck.lrnr
stck.lrnr.tuned <- setLearnerId(stck.lrnr.tuned, "stacked_learner")
ParamHelpers::getParamSet(stck.lrnr)

# make a parameter set to explore
hyperparams_to_tune <- makeParamSet(
  makeNumericParam("regr.randomForest.ntree", lower = 1, upper = 1e2),
  makeNumericParam("regr.glmnet.lambda", lower = 0, upper = 1e2)
)

# Make a Grid to Search On
ctrl <- makeTuneControlGrid()
# Choose a performance measure
optmeas <- rmse
# resampling approach with spatial CV considered
rdesc <-makeResampleDesc("SpRepCV", fold = 5, reps = 5)

ret <- tuneParams(learner = stck.lrnr, 
                  task = regr.task, 
                  resampling = rdesc, 
                  par.set = hyperparams_to_tune,
                  control = ctrl,
                  measures = optmeas, 
                  show.info = FALSE)

# apply tuning results
stck.lrnr.tuned <- setHyperPars(stck.lrnr, 
                                regr.randomForest.ntree = ret$x$regr.randomForest.ntree, 
                                regr.glmnet.lambda = ret$x$regr.glmnet.lambda
                                )
stck.lrnr.tuned <- setLearnerId(stck.lrnr.tuned, "stacked_learner_tuned")

# benchmark experiment
bench.lrnrs <-  list(stck.lrnr, stck.lrnr.tuned)
bmr <- benchmark(learners = bench.lrnrs, tasks = regr.task, resamplings = rdesc, measures = optmeas, 
                show.info = FALSE)











