library(tidyverse)
library(SuperLearner)
library(mlr)
library(mlrwrapSL)
source("R/00-IPTW_functions.R")


#............................................................
# Get Data and Setup
#............................................................
# simulate data from multivariate normal and make confounders
# Thanks, Lucy D'Agostino McGowan -- https://livefreeordichotomize.com/2019/01/17/understanding-propensity-score-weighting/
set.seed(48)
n <- 1000
X <- mvtnorm::rmvnorm(n,
                      mean = c(0.5, 1),
                      sigma = matrix(c(2, 1, 1, 1), ncol = 2))
dat.full <- tibble(
  x_1 = X[, 1],
  x_2 = X[, 2],
  treatment = as.numeric(- 0.5 + 0.25 * x_1 + 0.75 * x_2 + rnorm(n, 0, 1) > 0),
  outcome = 2 * treatment + rnorm(n, 0, 1))

glimpse(dat.full)
# drop outcome for now since we are fitting tx
dat <- dat.full %>% 
  dplyr::select(-c("outcome"))



#..............................
# Validation Sets etc. 
#.............................. 
# 50/50 split
nset <- 1:nrow(dat)
valset <- sort( sample(nset, nrow(dat)/2) )
valset.fct <-  ifelse(nset %in% valset, 1, 2) 
valset.list <- split(nset, factor(valset.fct))

#..............................
# Test, Train, Task objects
#.............................. 
x_train = dat[valset, c("x_1", "x_2")]
x_val = dat[-valset, c("x_1", "x_2")]
y_train = unlist( dat[valset, "treatment"] )
y_val = unlist( dat[-valset, "treatment"] )


#..............................
# SuperLearner
#.............................. 
learnerlib.SL <- c("SL.lm", "SL.glmnet", "SL.ksvm", "SL.ranger") # note, SL.glmnet has made this a lasso
sl = SuperLearner(Y = y_train, X = x_train,
                  SL.library = learnerlib.SL)

summary( sl$Z )
summary( sl$library.predict )

sl.results.Z <- sl$Z
sl.results.Zprime <- sl$library.predict
sl$coef

#..............................
# My Superlearner
#.............................. 
baselearners.names <- c("regr.lm", "regr.glmnet", "regr.ksvm", "regr.ranger")
baselearners <- lapply(baselearners.names, mlr::makeLearner)
baselearners[[2]] <- mlr::setHyperPars(baselearners[[2]] , alpha = 1) # make like SL above

task <- mlr::makeRegrTask(data = dat[valset, ], target = "treatment")
mySL <- mlrwrapSL::SL_crossval_risk_pred(learnerlib = baselearners, task = task, valset.list = list(c(1:250), c(250:500)))

mysl.results.Z <- mySL$Z
mysl.results.Zprime <- mySL$Zprime
mySL$cvrisk.coef

#..............................
# Look At Results
#.............................. 
baselearners.names <- gsub("regr.", "", baselearners.names)
sl.results.Z <- sl.results.Z %>% 
  tibble::as.tibble(.) %>% 
  magrittr::set_colnames(baselearners.names) %>% 
  dplyr::mutate(algorithm = "SuperLearner",
                level = "Z")

sl.results.Zprime <- sl.results.Zprime %>% 
  tibble::as.tibble(.) %>% 
  magrittr::set_colnames(baselearners.names) %>% 
  dplyr::mutate(algorithm = "SuperLearner",
                level = "Zprime")


mysl.results.Z <- mysl.results.Z %>% 
  tibble::as.tibble(.) %>% 
  magrittr::set_colnames(baselearners.names) %>% 
  dplyr::mutate(algorithm = "NFBSuperLearner",
                level = "Z")

mysl.results.Zprime <- mysl.results.Zprime %>% 
  tibble::as.tibble(.) %>% 
  magrittr::set_colnames(baselearners.names) %>% 
  dplyr::mutate(algorithm = "NFBSuperLearner",
                level = "Zprime")


SLresults <- rbind.data.frame(sl.results.Z, sl.results.Zprime, mysl.results.Z, mysl.results.Zprime)
SLresults %>% 
  dplyr::select(c("algorithm", "level", dplyr::everything())) %>% 
  tidyr::gather(., key = "baselearner", value = "pred", 3:ncol(.)) %>% 
  ggplot() + 
  geom_boxplot(aes(x=baselearner, y=pred, fill=algorithm)) + 
  facet_wrap(~level) + 
  ylab("Predictions") + ("Base Learner")











