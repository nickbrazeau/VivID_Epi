source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) 
library(rgeos)
library(rgdal)
library(geoR)
library(raster)
library(PrevMap)
set.seed(48)

#......................
# Import Data
#......................
gp.mod.framework <- readRDS("analyses/07-spatial_prediction/_rslurm_Prevmap_Long_Chain/params.RDS")

mcmc.ret <- list.files(path = "analyses/07-spatial_prediction/_rslurm_Prevmap_Long_Chain/", 
                       full.names = T, pattern = ".RDS")
mcmc.ret <- mcmc.ret[!grepl("params|f.", mcmc.ret)]

# 0-1 index here so don't have to worry about char 10 and 1 being next to eachother
gp.mod.framework$mcmc <- purrr::map(mcmc.ret, readRDS)

#...............................
# make prediction surfaces
#...............................

# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long","lat")

gp.mod.framework$grid.pred <- lapply(1:nrow(gp.mod.framework), function(x) return(grid.pred))

pred_PrevMap_bayes_wrapper <- function(mcmc, grid.pred){
  ret <- PrevMap::spatial.pred.binomial.Bayes(object = mcmc, grid.pred = grid.pred,
                                              type = "marginal", 
                                              scale.predictions = "prevalence",
                                              quantiles = NULL, 
                                              standard.error = T, 
                                              thresholds = NULL)
  return(ret)
  
}

####################################################################################
###########                           Rslurm                               #########
#####################################################################################

paramsdf <- gp.mod.framework[, c("mcmc", "grid.pred")]
setwd("analyses/07-spatial_prediction")
ntry <- nrow(mod.framework)
sjob <- rslurm::slurm_apply(f = pred_PrevMap_bayes_wrapper, 
                            params = paramsdf, 
                            jobname = 'Prevmap_predictions',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))

