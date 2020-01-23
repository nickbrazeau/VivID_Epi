#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get spatial predictions from 
# spatial models produced from prevmap
#----------------------------------------------------------------------------------------------------
source("R/00-functions_maps.R") 
source("R/00-functions_basic.R")
library(tidyverse)
library(raster)
library(PrevMap)
set.seed(48, "L'Ecuyer")

#......................
# Import Data
#......................
intercept <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-intercept.rds")
covar <- readRDS("analyses/07-spatial_prediction/prevmap_long_chains/longrun-prevmapbayes-mod.rds")
# to do add covars if needed
gp.mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                   mcmc = list(intercept, covar))


#...............................
# sample coordinates for prediction surface 
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.coords <- splancs::gridpts(poly, xs=0.05, ys=0.05)
colnames(grid.pred.coords) <- c("longnum","latnum")
grid.pred.coords.df <- as.data.frame(grid.pred.coords)

#...............................
# extract covariate information for prediction surface 
#...............................
# Raster surfaces for risk factors
riskvars = c("precip_mean_cont_scale_clst", 
             "cropprop_cont_scale_clst", "nightlightsmean_cont_scale_clst")
precipraster <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 

# take mean here since we converted to binary
cropraster <- readRDS("data/derived_data/vividepi_cropland_surface.rds")
cropraster <- raster::aggregate(cropraster, fact = 18, fun = mean)

# sum here under the assumption that night lights are additive (population density is additive)
nightligthraster <- raster::raster("data/derived_data/vividepi_nightlights_surface.grd")
nightligthraster <- raster::aggregate(nightligthraster, fact = 12, fun = sum)

# stack 
predcovars <- raster::stack(precipraster, cropraster, nightligthraster)
names(predcovars) <- riskvars

pred.df <- raster::extract(
  x = predcovars,
  y = sf::as_Spatial(
    sf::st_as_sf(grid.pred.coords.df, coords = c("longnum", "latnum"), 
                 crs = "+proj=utm +zone=34 +datum=WGS84 +units=m")),
  buffer = 6000,
  fun = mean,
  na.rm = T,
  sp = F)

pred.df <- as.data.frame(pred.df)
pred.df <- cbind.data.frame(grid.pred.coords.df, pred.df)

#...............................
# Bounds so we only do interpolation
# (i.e. coerce extrapolation values)
#...............................
# precip
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
precip <- dt %>% 
  dplyr::select("hv001", "precip_mean_cont_clst") %>% 
  dplyr::filter(!duplicated(.))
max.precip <- max(precip$precip_mean_cont_clst)
min.precip <- min(precip$precip_mean_cont_clst)

pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst > max.precip] <- max.precip
pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst < min.precip] <- min.precip


# bring in cropland
crop <- readRDS("data/derived_data/vividepi_cropland_propmeans.rds") 
max.crop <- max(crop$cropprop)
min.crop <- min(crop$cropprop)
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst > max.crop] <- max.crop
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst < min.crop] <- min.crop


# bring in nightlights
nightlists <- readRDS("data/derived_data/vividepi_night_clstmeans.rds") 
max.night <- max(nightlists$nightlightsmean)
min.night <- min(nightlists$nightlightsmean)
pred.df$nightlightsmean_cont_scale_clst[pred.df$nightlightsmean_cont_scale_clst > max.night] <- max.night
pred.df$nightlightsmean_cont_scale_clst[pred.df$nightlightsmean_cont_scale_clst < min.night] <- min.night


#...............................
# Scale
#...............................
pred.df <- apply(pred.df, 2, my.scale) 

#..............................................................
# Convert to raster for sampling
#..............................................................
covar.rstr.pred <- raster::rasterFromXYZ(pred.df, 
                                         res = c(0.05, 0.05),
                                         crs = "+proj=utm +zone=34 +datum=WGS84 +units=m")

# save this out for later
saveRDS(covar.rstr.pred, file = "data/derived_data/vividepi_spatial_covar_feature_engineer.rds")

# down sample
covar.rstr.pred.downsmpl <- raster::sampleRandom(covar.rstr.pred, 
                                                 size = 1e3,
                                                 na.rm = T,
                                                 xy = T,
                                                 sp = F)

colnames(covar.rstr.pred.downsmpl)[1:2] <- c("longnum", "latnum")

#...............................
# Setup Map Dataframe  
#...............................
# set up grid.pred
gp.mod.framework$grid.pred <- list(covar.rstr.pred.downsmpl[,c("longnum", "latnum")],
                                   covar.rstr.pred.downsmpl[,c("longnum", "latnum")])

# set up predictors
gp.mod.framework$predictors <- list(NULL, covar.rstr.pred.downsmpl[,c("precip_mean_cont_scale_clst", 
                                                                      "cropprop_cont_scale_clst",
                                                                      "nightlightsmean_cont_scale_clst")])



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
                            jobname = 'Prevmap_predictions_large',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = "500g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))




