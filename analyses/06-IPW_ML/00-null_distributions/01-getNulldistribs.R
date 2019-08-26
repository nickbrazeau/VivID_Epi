#----------------------------------------------------------------------------------
# Purpose of this script is to make the null distributions of the covariate maps
#-------------------------- --------------------------------------------------------

source("R/00-functions_epi.R")
source("R/00-make_null_IPTW_distribs_brownian.R")
source("R/00-functions_Ensemble_Wrapper.R")
set.seed(48, "L'Ecuyer")
library(tidyverse)
library(mlr)

dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL



#............................................
# read in treatment map file which 
# are our covariates/txs to be estimated
#............................................
# read in treatments
txs <- readRDS("model_datamaps/IPTW_treatments.RDS") %>% 
  dplyr::rename(positive = positivefactor) %>% 
  dplyr::mutate(type = ifelse(grepl("cont", column_name), "continuous",
                              ifelse(grepl("fctb", column_name), "binary", NA)))


#........................
# manipulate data
#........................
# subset to treatments, outcome, weights and coords
varstoinclude <- c("pv18s" , "pfldh", "hv005_wi", txs$column_name,
                   "alt_dem_cont_scale_clst", "hab1_cont_scale", "hv104_fctb", "wtrdist_cont_log_scale_clst", # need to add in covariates that don't have confounding ancestors but are needed elsewhere
                   "hiv03_fctb", # no longer considered risk factor bc too few observations
                   "longnum", "latnum")
dt.ml <- dt %>% 
  dplyr::select(varstoinclude)

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

txs <- txs %>% 
  dplyr::select(c("target", "data", "adj_set"))


#................................
# Build dataframe for rSLURM
#................................

nulliters <- 1e3
paramsdf <- lapply(1:nulliters, function(x) return(txs)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(target) %>% 
  dplyr::rename(covars = adj_set)



# for slurm on LL
setwd("analyses/06-IPW_ML/00-null_distributions/")
ntry <- nulliters
sjob <- rslurm::slurm_apply(f = make.null.distribution.energy, 
                            params = paramsdf, 
                            jobname = 'nulldist_vivid_covar',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 64000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 128),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))




