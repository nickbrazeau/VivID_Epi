source("R/00-simple_Ensemble_Wrapper.R")
library(tidyverse)
library(mlr)

dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
#........................
# manipulate tx map
#........................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(target = column_name,
                positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", target), "continuous",
                              ifelse(grepl("fctb", target), "binary", NA)))
  
#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
dt.ml <- dt %>% 
  dplyr::select(c("pv18s" , "pfldh", "hv005_wi", txs$target, 
                "longnum", "latnum"))

# subset to complete cases
dt.ml <- dt.ml %>% 
  dplyr::filter(complete.cases(.)) %>% 
  data.frame(.)


#........................
# Subset and Store Dataframes for Tasks
#........................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml %>% 
                            dplyr::select(c(x,y))
                          return(ret)
                        })



#-------------------------------------------------
# Run the Models
#-------------------------------------------------

# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type")], 
                         .f = make_class_task)

# now look and correct class imbalance
txs$task <- purrr::pmap(txs[,c("task", "type")], 
                         .f = find_Class_Imbalance,
                        classimb_tol = 0.6,
                        smotenn = 5)


# now make the ensemble learner
txs$learner <- purrr::map(txs$type, make_simple_Stack, 
                          learners = baselearners.list)

# make full model
txs$fullmodel <- purrr::map2(.x = txs$learner, .y = txs$task,
                             train)


# get predictions
txs$preds <- purrr::map2(.x = txs$fullmodel, .y = txs$task,
                             predict)


saveRDS(object = txs, file = "results/ensemble_predictive_probs.RDS")





