#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get spatial predictions from 
# spatial models produced from prevmap
#----------------------------------------------------------------------------------------------------
source("R/00-functions_maps.R") 
library(tidyverse)
library(raster)
library(PrevMap)
set.seed(48)

#......................
# Import Data
#......................
intercept <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-intercept.rds")
covar <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-mod.rds")
# to do add covars if needed
gp.mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                   mcmc = list(intercept, covar))


#...............................
# make prediction surfaces for intercept
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.intercept <- splancs::gridpts(poly, xs=0.05, ys=0.05)
colnames(grid.pred.intercept) <- c("longnum","latnum")


# Raster surfaces for risk factors
riskvars = c("precip_mean_cont_scale_clst", "temp_mean_cont_scale_clst", 
             "cropprop_cont_scale_clst", "nightlightsmean_cont_scale_clst")
precipraster <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 
tempraster <- readRDS("data/derived_data/vividepi_temperature_study_period_effsurface.rds")
cropraster <- readRDS("data/derived_data/vividepi_cropland_surface.rds")
nightlisthraster <- raster::raster("data/derived_data/vividepi_nightlights_surface.grd")

# reproject so all on same scale
# precipraster
# tempreaster
cropraster.repr <- raster::projectRaster(cropraster, precipraster)
nightlisthraster.repr <- raster::projectRaster(nightlisthraster, precipraster)


# manipulate crop raster
xy <- raster::coordinates(cropraster.repr)
xy.list <- split(xy, 1:nrow(xy))
xy.list <- lapply(xy.list, function(x){
  ret <- sf::st_sfc(sf::st_point(x), 
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  return(ret)
})


cropraster.smooth.values <- sapply(xy.list, function(sppoint){
  ret <- raster::extract(x = cropraster,
                  y = sf::as_Spatial(sppoint),
                  buffer = 6000,
                  fun = mean,
                  na.rm = T,
                  sp = F)
  return(ret)
})

# overlay new values
cropraster.smooth <- cropraster.repr
values(cropraster.smooth) <- cropraster.smooth.values

# get predictive df
pred.df <- data.frame(longnum = raster::coordinates(precipraster)[,1],
                      latnum = raster::coordinates(precipraster)[,2],
                      precip_mean_cont_scale_clst = raster::values(precipraster),
                      temp_mean_cont_scale_clst = raster::values(tempraster),
                      cropprop_cont_scale_clst = raster::values(cropraster.smooth),
                      nightlightsmean_cont_scale_clst = raster::values(nightlisthraster.repr)
) %>% 
  dplyr::filter(!is.na(precip_mean_cont_scale_clst),
                !is.na(temp_mean_cont_scale_clst),
                !is.na(cropprop_cont_scale_clst),
                !is.na(nightlightsmean_cont_scale_clst)) 

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
                                              quantiles = NULL, 
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
                            slurm_options = list(mem = "500g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))
