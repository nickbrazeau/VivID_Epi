source("R/00-functions_epi.R")
source("R/00-make_null_IPTW_distribs_brownian.R")
source("R/00-my_IPTW_performance_measure_energy.R")
set.seed(48, "L'Ecuyer")
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
# Setup Models
#-------------------------------------------------

# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type")], 
                        .f = make_class_task)

# now look and correct class imbalance
txs$task <- purrr::pmap(txs[,c("task", "type")], 
                        .f = find_Class_Imbalance,
                        classimb_tol = 0.6,
                        smotenn = 5)
txs %>% 
  dplyr::select(c("target", "task", "adj_set", "data"))

nulliters <- 1e3
txs$nulldist <- purrr::pmap(txs, 
                            function(target, task, adj_set){
                              
                              data <- mlr::getTaskData(task)
                              
                              # R for loop
                              ret <- parallel::mclapply(1:nulliters, function(x){ 
                                return( make.null.distribution.energy(data = data, covars = adj_set, target = target) )
                                
                              }) %>% unlist(.)
                              
                            })

saveRDS(txs, file = "analyses/06-IPW_ML/nulldist.RDS")

