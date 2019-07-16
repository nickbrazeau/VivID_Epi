#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make my own performance measure to estimate how well
# our IPTWeights are balancing the baseline covariates
# This is inspired by/is a slight extension of the function/methods presented in 
# Y. Zhu et al. "A Boosting Algorithm ..." 2015
#----------------------------------------------------------------------------------------------------
source("R/00-functions_epi.R")
source("R/00-make_null_IPTW_distribs.R")

#..............................................................
# Binary Treatment Case
#..............................................................
my.covarbal.binary.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame for Dij calculations
  dat <- mlr::getTaskData(task)
  n <- nrow(dat)
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  # get inverse probability weights
  wi <- get_iptw_prob(task = task, preds = pred, type = "binary")
  
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



# Generate the Measure object for binary tx
my.covarbal.binary = makeMeasure(
  id = "my.covarbal.binary", 
  name = "Baseline covariate balance estimator for binary tx case",
  properties = c("classif"),
  extra.args = list(nulldist = NA),
  minimize = TRUE, 
  best = -Inf, worst = Inf,
  fun = my.covarbal.binary.fun
)



#..............................................................
# Continuous Tx Case
#..............................................................
my.covarbal.continuous.fun = function(task, model, pred, feats, nulldist) {
  
  # pull in pieces to make data frame for Dij calculations
  dat <- mlr::getTaskData(task)
  n <- nrow(dat)
  target <- mlr::getTaskTargetNames(task)
  covars <- mlr::getTaskFeatureNames(task)
  
  # get inverse probability weights
  wi <- get_iptw_prob(task = task, preds = pred, type = "continuous")
  
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
  dat.map$cor <- purrr::map(dat.map$widata, my.pearson, target = target)
  # average dj
  avgdj <- mean( unlist(dat.map$cor) )
  
  # now calculate Z score
  Z <- ( avgdj - mean(nulldist) ) / sd(nulldist)
  
  return(Z)
  }



# Generate the Measure object for continuous tx
my.covarbal.continuous = makeMeasure(
  id = "my.covarbal.binary", 
  name = "Baseline covariate balance estimator for continuous tx",
  properties = c("regr"),
  extra.args = list(nulldist = NA),
  minimize = TRUE, 
  best = -Inf, worst = Inf,
  fun = my.covarbal.continuous.fun
)


































