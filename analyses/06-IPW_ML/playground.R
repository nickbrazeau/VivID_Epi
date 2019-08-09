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


source("R/00-make_null_IPTW_distribs.R")
library(tidyverse)
set.seed(928)
n <- 1000
X <- mvtnorm::rmvnorm(n,
                      mean = c(0.5, 1),
                      sigma = matrix(c(2, 1, 1, 1), ncol = 2)
)

dat <- tibble(
  x_1 = X[, 1],
  x_2 = X[, 2],
  treatment = as.numeric(- 0.5 + 0.25 * x_1 + 0.75 * x_2 + rnorm(n, 0, 1) > 0),
  outcome = 2 * treatment + rnorm(n, 0, 1)
) %>% 
  dplyr::mutate(treatment = factor(treatment))
glimpse(dat)


null.binaryTx <- sapply(1:1e4, function(x){
  make.null.distribution.binaryTx(target = "treatment", covars = c("x_1", "x_2"), data = dat)})

hist(null.binaryTx)


target = "treatment"
covars = c("x_1", "x_2")

# task to f it
task = makeClassifTask(data = dat[c(target, covars)], target = "treatment", positive = "1")
lrn <- makeLearner("classif.glmnet", predict.type = "prob")
fullmodel <- mlr::train(learner = lrn,
                        task = task)
predictions <- predict(fullmodel, task = task)




# Define a function that calculates the misclassification rate
my.covarbal.binary.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame
  dat <- mlr::getTaskData(task)
  n <- nrow(dat)
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  wi <- mlr::getPredictionProbabilities(pred)
  
  # coerce back to vector
  nulldist <- unname( unlist(nulldist) )
  
  if(!is.atomic(nulldist)){
    stop("There is an issue with coercion from mlr performance to 
         a vector that we can call")
  }
  
  
  # pull apart data to find Djs
  dat.map <- dat %>% 
    dplyr::select(c(target, covars)) %>% 
    tidyr::gather(., key = "covar", value = "val", covars) %>% 
    dplyr::group_by(covar) %>% 
    tidyr::nest()
  
  dat.map$widata <- purrr::map(dat.map$data, function(x){
    x <- x[sample(1:nrow(x), size = n, prob = wi, replace = T), ]
    return(x)
  })
  dat.map$smd <- purrr::map(dat.map$widata, smd, target = target)
  # average dj
  avgdj <- mean( unlist(dat.map$smd) )
  
  # now calculate Z score
  Z <- ( avgdj - mean(nulldist) ) / sd(nulldist)
  
  return(Z)
  }



# Generate the Measure object
my.covarbal.binary = makeMeasure(
  id = "my.covarbal.binary", 
  name = "Baseline covariate balance estimator",
  properties = c("classif"),
  extra.args = list(nulldist = NA),
  minimize = TRUE, 
  best = -Inf, worst = Inf,
  fun = my.covarbal.binary.fun
)

my.covarbal.binary <- mlr::setMeasurePars(my.covarbal.binary, 
                                          par.vals = list(nulldist = null.binaryTx))


performance(pred = predictions, measures = my.covarbal.binary, 
            task = task)




ps = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 1))

rdesc = makeResampleDesc("CV", iters = 3L)
# # 3) Here we use Random Search (with 10 Iterations) to find the optimal hyperparameter
ctrl =  makeTuneControlRandom(maxit = 10)
# # 4) now use the learner on the training Task with the 3-fold CV to optimize your set of parameters and evaluate it with SQWK
res = tuneParams(lrn, 
                 task = task, 
                 resampling = rdesc, 
                 par.set = ps, 
                 control = ctrl, 
                 measures = my.covarbal.binary
)











ps = makeParamSet(
  makeNumericParam("alpha", lower = 0, upper = 1))

rdesc = makeResampleDesc("CV", iters = 3L)
# # 3) Here we use Random Search (with 10 Iterations) to find the optimal hyperparameter
ctrl =  makeTuneControlRandom(maxit = 10)
# # 4) now use the learner on the training Task with the 3-fold CV to optimize your set of parameters and evaluate it with SQWK
res = tuneParams(lrn, 
                 task = task, 
                 resampling = rdesc, 
                 par.set = ps, 
                 control = ctrl, 
                 measures = my.covarbal.binary
)












Housing2 %>% 
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
stck.lrnr <- setLearnerId(stck.lrnr, "stacked_learner")
ParamHelpers::getParamSet(stck.lrnr)

# make a parameter set to explore
hyperparams_to_tune <- makeParamSet(
  makeNumericParam("regr.glmnet.alpha", lower = 0, upper = 1),
  makeNumericParam("regr.kknn.k", lower = 1, upper =5 ),
  makeNumericParam("regr.ksvm.C", lower = 1, upper = 3),
  makeDiscreteParam("regr.ksvm.kernel", values = c("rbfdot", "vanilladot")), # values = c("rbfdot", "polydot", "vanilladot", "tanhdot", "laplacedot", "besseldot", "anovadot", "splinedot", "stringdot")),
  makeNumericParam("regr.randomForest.mtry", lower = 1, upper = 3 )
)

# L1/L2 Regularization (glmnet), alpha: The elasticnet mixing parameter, with 0≤α≤ 1. The penalty is defined as (1-α)/2||β||_2^2+α||β||_1 alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
# K-Nearest Neighbors, k: Number of neighbors considered.
# Single Vector Machine, C: cost of constraints violation (default: 1) this is the `C'-constant of the regularization term in the Lagrange formulation.
# Single Vector Machine, kernel: The kernel function used in training and predicting. This parameter can be set to any function, of class kernel, which computes the inner product in feature space between two vector arguments (see kernels). kernlab provides the most popular kernel functions which can be used by setting the kernel parameter to the following strings: 
# Random Forest, Mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 




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
                  show.info = T)

# apply tuning results
stck.lrnr.tuned <- setHyperPars(stck.lrnr, 
                                regr.glmnet.alpha = ret$x$regr.glmnet.alpha,
                                regr.kknn.k = ret$x$regr.kknn.k,
                                regr.ksvm.C = ret$x$regr.ksvm.C,
                                regr.ksvm.kernel = ret$x$regr.ksvm.kernel,
                                regr.randomForest.mtry = ret$x$regr.randomForest.mtry
                                )
stck.lrnr.tuned <- setLearnerId(stck.lrnr.tuned, "stacked_learner_tuned")

# benchmark experiment
bench.lrnrs <-  list(stck.lrnr, stck.lrnr.tuned)
bmr <- benchmark(learners = bench.lrnrs, 
                 tasks = regr.task, 
                 resamplings = rdesc, 
                 measures = optmeas, 
                show.info = T)

getBMRAggrPerformances(bmr)
plotBMRBoxplots(bmr)

pred <- train(stck.lrnr, regr.task)








