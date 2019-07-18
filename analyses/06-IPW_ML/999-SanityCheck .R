#----------------------------------------------------------------------------------------------------
# purpose of this script is to compare and ensemble model 
# from SuperLearner and MLR Ensemble SuperLearner
# For the SuperLearner Package, will follow their Vignette (https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html)
# then do the same thing in MLR
#----------------------------------------------------------------------------------------------------
library(tidyverse)

#................................................
# Data Setup (following Superlearner vignette)
#................................................
library(MASS)
data(Boston)
?Boston
glimpse(Boston)

# Extract our outcome variable from the dataframe.
outcome = Boston$medv

# Create a dataframe to contain our explanatory variables.
data = subset(Boston, select = -medv)

#................................................
# Train/Subset Setup for Superlearners
#................................................
# Set a seed for reproducibility in this random sampling.
set.seed(1)

# Reduce to a dataset of 150 observations to speed up model fitting.
train_obs = sample(nrow(data), 150)

# X is our training sample.
X_train = data[train_obs, ]

# Create a holdout set for evaluating model performance.
# Note: cross-validation is even better than a single holdout sample.
X_holdout = data[-train_obs, ]

# Create a binary outcome variable: towns in which median home value is > 22,000.
outcome_bin = as.numeric(outcome > 22)

Y_train = outcome_bin[train_obs]
Y_holdout = outcome_bin[-train_obs]



#-----------------------------------------------------------------------
# SUPERLEARNER PACKAGE
#-----------------------------------------------------------------------
#remotes::install_github("ecpolley/SuperLearner")
library(SuperLearner)

sl.trainmodel <- SuperLearner(Y = Y_train, X = X_train, family = binomial(),
                              SL.library = c("SL.glmnet", "SL.randomForest"))
sl.trainmodel

sl.pred = predict(sl.trainmodel, X_holdout, onlySL = T)

# Check the structure of this prediction object.
str(sl.pred)
# Review the columns of $library.predict.
summary(sl.pred$library.predict)
# Look at actual predictions
summary(sl.pred$pred)


#-----------------------------------------------------------------------
# MLR PACKAGE
#-----------------------------------------------------------------------
#remotes::install_github("mlr-org/mlr")
library(mlr)

# setup mlr task
mlr.data <- cbind.data.frame(Y = outcome_bin, data) %>% 
  dplyr::mutate(Y = factor(Y))
test_obs <- seq(1:nrow(data))[ !seq(1:nrow(data)) %in% train_obs ]

# setup stacked learner, e.g. specificy libraries
mlr.task <- makeClassifTask(id = "mlrtest", data = mlr.data, target = "Y", positive = "1")
baselearners <- c("classif.glmnet", "classif.randomForest")
baselearners <- lapply(baselearners, makeLearner, predict.type = "prob")
mlr.stck.lrnr <- makeStackedLearner(base.learners = baselearners,
                                    predict.type = "prob", 
                                    method = "hill.climb")
                                  
                                  #  super.learner = "classif.lda",
                                  #  method = "stack.cv",
                                  #  resampling = makeResampleDesc("RepCV", fold = 5, reps = 5))

mlr.trainmodel <- mlr::train(learner = mlr.stck.lrnr, task = mlr.task, subset = train_obs)
mlr.pred = predict(mlr.trainmodel, task = mlr.task, subset = test_obs)
mlr.pred



#................................................
# Compare Performance
#................................................
# Review AUC - Area Under Curve
sl.pred_rocr = ROCR::prediction(sl.pred$pred, Y_holdout)
sl.auc = ROCR::performance(sl.pred_rocr, measure = "auc", x.measure = "cutoff")@y.values[[1]]
sl.auc

mlr::performance(mlr.pred, measures = list(mlr::auc))


summary(unlist( sl.pred_rocr@predictions ))
summary(mlr.pred$data$prob.1)
plot(unlist( sl.pred_rocr@predictions ) ~ mlr.pred$data$prob.1, main = "SL Predictive Prob vs. MLR Predictive Prob")

# it apepars that the Hill-Climbing Algorithm Approximates the 
# nnls meta function in the SuperLearner R package





