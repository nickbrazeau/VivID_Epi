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












