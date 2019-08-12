#----------------------------------------------------------------------------------------------------
# Purpose of this script is to apply Tuning Parameters and Make Final Models 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(mlr)
source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")

# after figuring out who won
# this will get hyperpar mlr::getHyperPars(learner.hp)



dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
#........................
# manipulate tx map
#........................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA)))
# note we want the targets to be in their original form and no the scaled form (to account for variance in the outcome, which is now our tx level)
txs <- txs %>% 
  dplyr::mutate(target = column_name,
                target = gsub("_scale", "", column_name))


#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
varstoinclude <- c("pv18s" , "pfldh", "hv005_wi", txs$target, txs$column_name,
                   "alt_dem_cont_scale_clst", "hab1_cont_scale", "hv104_fctb", "wtrdist_cont_scale_clst", # need to add in covariates that don't have confounding ancestors but are needed elsewhere
                   "urbanscore_cont_scale_clst", #urbanicity
                   "longnum", "latnum")

dt.ml <- dt %>% 
  dplyr::select(varstoinclude)

# subset to complete cases
dt.ml.cc <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum")) %>% 
  data.frame(.)

dt.ml.coords <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  data.frame(.)


#........................
# Subset and Store Dataframes for Tasks
#........................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml.cc %>% 
                            dplyr::select(c(x, y))
                          return(ret)
                        })

#........................
# Pull Out Coords
#........................
txs$coordinates <- lapply(1:nrow(txs), function(x) return(dt.ml.coords))


#--------------------------------------
# Setup tasks & base learners
#--------------------------------------
# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type", "coordinates")], 
                        .f = make_class_task)


txs$learner <- purrr::map(txs$task, make_avg_Stack, 
                          learners = baselearners.list)


#--------------------------------------
# Train & Predictions
#--------------------------------------
txs$subset <- purrr::map(txs$task, function(x){return( sample(mlr::getTaskSize(x), size = n/3) )})

# Train the learner
txs$train <-  purrr::pmap(txs[,c("learner", "task", "subset")], mlr::train)

# get predictions
txs$preds <- purrr::pmap(txs[,c("task", "train")], function(task, train){ ret <- predict(train, task); return(ret) })


# make inverse probability weights
txs$iptw <- pmap(txs[,c("task", "preds")], get_iptw_prob)


# analyze weights
lapply(txs$iptw, summary) %>% 
  do.call("rbind.data.frame", .) %>% 
  dplyr::mutate_if(is.numeric, round, 2) %>% 
  magrittr::set_colnames(c("min", "1stquart", "median", "mean", "3rdqart", "max"))



# trainLearner


temp <- txs[2,]
View( cbind( temp$preds[[1]]$data, temp$iptw[[1]]))


