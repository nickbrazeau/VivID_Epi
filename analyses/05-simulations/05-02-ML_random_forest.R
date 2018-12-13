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
source("analyses/05-simulations/05-00-Create_Sim_Data.R")
