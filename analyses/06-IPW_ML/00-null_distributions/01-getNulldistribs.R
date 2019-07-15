#----------------------------------------------------------------------------------
# Purpose of this script is to make the null distributions of the covariate maps
#----------------------------------------------------------------------------------

source("R/00-functions_epi.R")
source("R/00-make_null_IPTW_distribs_brownian.R")
source("R/00-simple_Ensemble_Wrapper.R")
set.seed(48, "L'Ecuyer")
library(tidyverse)
library(mlr)

dt <- readRDS("data/derived_data/vividepi_recode.rds")
#sf::st_geometry(dt) <- NULL



#............................................
# read in treatment map file which 
# are our covariates/txs to be estimated
#............................................
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


#.........................................
# Subset and Store Dataframes for Tasks
#.........................................
txs$data <- purrr::map2(.x = txs$target, .y = txs$adj_set, 
                        .f = function(x, y){
                          ret <- dt.ml %>% 
                            dplyr::select(c(x,y))
                          return(ret)
                        })



#........................
# Setup Models
#........................

# first make the tasks
txs$task <- purrr::pmap(txs[,c("data", "target", "positive", "type")], 
                        .f = make_class_task)

# now look and correct class imbalance
txs$task <- purrr::pmap(txs[,c("task", "type")], 
                        .f = find_Class_Imbalance,
                        classimb_tol = 0.6,
                        smotenn = 5)
txs <- txs %>% 
  dplyr::select(c("target", "task", "adj_set"))


#................................
# Build dataframe for rSLURM
#................................

nulliters <- 1e3
paramsdf <- lapply(1:nulliters, function(x) return(txs)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(target)


slurm_get_nulldist <-  function(target, task, adj_set){
  
  data <- mlr::getTaskData(task)
  ret <- make.null.distribution.energy(data = data, covars = adj_set, target = target)
  return(ret)
  }
  

# for slurm on LL
setwd("analyses/06-IPW_ML/00-null_distributions/")
ntry <- 128
sjob <- rslurm::slurm_apply(f = slurm_get_nulldist, 
                            params = paramsdf, 
                            jobname = 'nulldist_vivid_covar',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry - 1, 
                                                                 16),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))




