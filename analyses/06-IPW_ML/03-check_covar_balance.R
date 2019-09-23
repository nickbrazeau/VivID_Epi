#----------------------------------------------------------------------------------------------------
# Purpose of this script is to check for 
# covariate balance after applying the weights
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(energy)



# pull covars twice for weights
# can just use the param table to get adj set and tx
# can make boxplots of the dij values 
# before and after you apply the weights
# using get task data