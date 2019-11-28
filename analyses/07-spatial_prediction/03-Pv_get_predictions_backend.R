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
# sample coordinates for prediction surface 
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.coords <- splancs::gridpts(poly, xs=0.05, ys=0.05)
colnames(grid.pred.coords) <- c("longnum","latnum")

#...............................
# downsample for computational burden
#...............................
coords_smpl <- sample(1:nrow(grid.pred.coords), 
                      size = 16000, replace = F)
coords_smpl <- sort(coords_smpl)
grid.pred.coords <- grid.pred.coords[coords_smpl, ]
grid.pred.coords.df <- as.data.frame(grid.pred.coords)
#...............................
# extract covariate information for prediction surface 
#...............................
# Raster surfaces for risk factors
riskvars = c("precip_mean_cont_scale_clst", "temp_mean_cont_scale_clst", 
             "cropprop_cont_scale_clst", "nightlightsmean_cont_scale_clst")
precipraster <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 
tempraster <- readRDS("data/derived_data/vividepi_temperature_study_period_effsurface.rds")
cropraster <- readRDS("data/derived_data/vividepi_cropland_surface.rds")
nightlisthraster <- raster::raster("data/derived_data/vividepi_nightlights_surface.grd")


extract_rstr_values <- function(rstr, coords){
  ret <- raster::extract(x = rstr,
                         y = sf::as_Spatial(sf::st_as_sf(coords, coords = c("longnum", "latnum"), 
                                                         crs = "+proj=utm +zone=34 +datum=WGS84 +units=m")),
                         buffer = 6000,
                         fun = mean,
                         na.rm = T,
                         sp = F)
  return(ret)
  
}
# note, every raster needs to be transformed 
rstrmeans <- lapply(list(precipraster, tempraster, cropraster, nightlisthraster), 
                    extract_rstr_values, coords = grid.pred.coords.df)

names(rstrmeans) <- riskvars
values(rstrmeans[["precip_mean_cont_scale_clst"]]) <- my.scale(values(rstrmeans[["precip_mean_cont_scale_clst"]]))
values(rstrmeans[["temp_mean_cont_scale_clst"]]) <- my.scale(values(rstrmeans[["temp_mean_cont_scale_clst"]]))
values(rstrmeans[["nightlightsmean_cont_scale_clst"]]) <- my.scale(values(rstrmeans[["nightlightsmean_cont_scale_clst"]]))
values(rstrmeans[["cropprop_cont_scale_clst"]]) <- my.scale(logit(values(rstrmeans[["cropprop_cont_scale_clst"]]), tol = 1e-3))


# get predictive df
pred.df <- data.frame(longnum = grid.pred.coords.df[,"longnum"],
                      latnum = grid.pred.coords.df[, "latnum"],
                      precip_mean_cont_scale_clst = values(rstrmeans[["precip_mean_cont_scale_clst"]]),
                      temp_mean_cont_scale_clst = values(rstrmeans[["temp_mean_cont_scale_clst"]]),
                      cropprop_cont_scale_clst = values(rstrmeans[["cropprop_cont_scale_clst"]]),
                      nightlightsmean_cont_scale_clst = values(rstrmeans[["nightlightsmean_cont_scale_clst"]])
) %>% 
  dplyr::filter(!is.na(precip_mean_cont_scale_clst),
                !is.na(temp_mean_cont_scale_clst),
                !is.na(cropprop_cont_scale_clst),
                !is.na(nightlightsmean_cont_scale_clst)) 



# set up grid.pred
gp.mod.framework$grid.pred <- list(grid.pred.coords, grid.pred.coords)

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
