# https://stackoverflow.com/questions/25715502/using-randomforest-package-in-r-how-to-get-probabilities-from-classification-mo
# https://shirinsplayground.netlify.com/2018/10/ml_basics_rf/
# https://www.r-bloggers.com/machine-learning-basics-gradient-boosting-xgboost/
# RF uses bagging (Bootstrap aggregation)
#   RF uses random subset of data (bagging) as well as at each node it performs feature bagging and choses a random subset of features to split on
# Gradient boosting machines are also made up of decision trees and use boosting (weighted sampling)
#   The idea behind the weights is that hard cases should come up more in learning/training sets -- by increasing the sampling of incorrectly predicted observations
# Goal is to have individual trees not be correlated (as to avoid bias) -- which is where random sampling comes into play

# Running random forests
# https://shirinsplayground.netlify.com/2018/06/intro_to_ml_workshop_heidelberg/
# 
# 

#......................
# Import Data and Dependencies
#......................
library(mlr)
load("~/Documents/GitHub/VivID_Epi/analyses/888-simulations/simdata/simdat_basic.rds")
expdat <- dat %>% 
  mutate(outcome = factor(outcome, levels = c(0,1), labels = c("neg", "pos")),
         treatment = factor(treatment, 
                            levels=c(0,1), labels = c("neg", "pos")),
         hv001 = factor(hv001)
  ) %>% 
  select(c("x_1", "x_2", "treatment")) %>% 
  as.data.frame(.)

xtabs(~expdat$treatment)
#......................
# Make the task (i.e. dataset) and the learner
#......................

classif.task = makeClassifTask(id = "EN_dat", 
                               data = expdat, 
                               target = "treatment")


#Classification tree, set it up for predicting probabilities
classif.rf = makeLearner("classif.randomForest", 
                         predict.type = "prob", 
                         fix.factors.prediction = TRUE)

getHyperPars(classif.rf)
getParamSet(classif.rf)

#......................
# Setup Train
#......................
n = getTaskSize(classif.task)
# Use 2/3 of the observations for training and 1/3 for testing
train.set = sample(n, size = 2*n/3)
test.set = seq(1,n)[ ! seq(1, n) %in% train.set ]
mod <- train(classif.rf, classif.task, subset = train.set)
task.pred <- predict(mod, task = classif.task, subset = test.set)
task.pred
calculateConfusionMatrix(task.pred)
# plotLearnerPrediction(classif.rf, task = classif.task)
performance(task.pred, measures = list(fpr, fnr, mmce))
d = generateThreshVsPerfData(task.pred, measures = list(fpr, fnr, mmce))
plotThreshVsPerf(d)
calculateROCMeasures(task.pred)




#......................
# IPW Weights
#......................
# dat$ps <- predict(mod, task = classif.task)$data$prob.pos
pred <- predict(mod, task = classif.task)$data
dat$ps <- ifelse(pred$response == "pos", pred$prob.pos, pred$prob.neg)
ggplot(data=dat,aes(x=ps, 
                    group=factor(treatment), fill=factor(treatment))) +
  geom_histogram(aes(y=..density..),alpha = 0.75,binwidth=0.02,position = position_dodge(width=0.01))+
  theme_classic()+
  xlab("Predicted probability of Ai")+
  labs(fill = "Observed")
# create weights
p_exposure <- sum(dat$treatment == "pos") / nrow(dat)
dat$iptw_s <-  1/dat$ps

weighted_df <- survey::svydesign(~0, weights = dat$iptw_s, data=dat)
survey::svymean(~treatment, weighted_df, na.rm=T)
tableone::svyCreateTableOne(vars = c("x_1", "x_2"), 
                  strata = "treatment", test = F, 
                  data = weighted_df)

dat$x <- 1:nrow(dat)
broom::tidy(geepack::geeglm(data = dat, outcome ~ treatment, 
                   family=binomial("logit"), 
                   weight=iptw_s, id = x))

