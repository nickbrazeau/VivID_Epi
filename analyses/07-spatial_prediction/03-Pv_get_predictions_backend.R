#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get spatial predictions from 
# spatial models produced from prevmap
#----------------------------------------------------------------------------------------------------
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
intercept <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-intercept.rds")
covar <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-mod.rds")
gp.mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                   mcmc = list(intercept, covar))


#...............................
# make prediction surfaces for intercept
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.intercept <- splancs::gridpts(poly, xs=0.5, ys=0.5)
colnames(grid.pred.intercept) <- c("longnum","latnum")


# note, because we have a raster surface of precipitation, we can 
# predict at all (or a random sample of) raster points

predictors <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 
# make this into a manageable size
predictors <- raster::sampleRandom(predictors, size = 1e4,
                                   asRaster = T)
pred.df <- data.frame(longnum = raster::coordinates(predictors)[,1],
                      latnum = raster::coordinates(predictors)[,2],
                      precip_mean_cont_scale_clst = raster::values(predictors)) %>% 
  dplyr::filter(!is.na(precip_mean_cont_scale_clst)) # this removes NAs introducted above
grid.pred.covars <- pred.df[,c("longnum", "latnum")]


# set up grid.pred
gp.mod.framework$grid.pred <- list(grid.pred.intercept, grid.pred.covars)

# set up predictors
gp.mod.framework$predictors <- list(NULL, pred.df)



# make wrapper
pred_PrevMap_bayes_wrapper <- function(mcmc, grid.pred, predictors){
  ret <- PrevMap::spatial.pred.binomial.Bayes(object = mcmc, 
                                              grid.pred = grid.pred,
                                              predictors = predictors, 
                                              type = "marginal", 
                                              scale.predictions = "prevalence",
                                              quantiles = c(0.025, 0.975), 
                                              standard.error = T, 
                                              thresholds = NULL)
  return(ret)
  
}

####################################################################################
###########                           Preds                                #########
#####################################################################################
paramsdf <- gp.mod.framework[, c("mcmc", "grid.pred", "predictors")]
setwd("analyses/07-spatial_prediction")
ntry <- nrow(gp.mod.framework)
sjob <- rslurm::slurm_apply(f = pred_PrevMap_bayes_wrapper, 
                            params = paramsdf, 
                            jobname = 'Prevmap_predictions',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 64000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))
