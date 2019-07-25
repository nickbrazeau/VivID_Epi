#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create a spatial prediction model
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) #wrap the survey package in dplyr syntax
library(PrevMap)

#......................
# Import Data
#......................
# Summarized by Cluster
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
mp <- mp %>% 
  dplyr::filter(maplvl == "hv001") 


#......................
# Subset to Pv
#......................
pvclust.weighted <- mp %>% 
  dplyr::filter(plsmdmspec == "pv18s" & maplvl == "adm1name") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pvclust.weighted <- sf::st_as_sf(pvclust.weighted)
# need to keep integers, so will round
pvprov.weighted <- pvprov.weighted %>% 
  dplyr::mutate_if(is.numeric, round, 0)

#----------------------------------------------------------------------------------------------------
# Smoothed Guassian Maps
#----------------------------------------------------------------------------------------------------





coords <- as.formula(paste0("~", long_var, "+", lat_var))
ret.fit <- PrevMap::linear.model.MLE(formula=eq, coords=coords, 
                                     data=data, start.cov.pars=start.cov.pars, 
                                     kappa=kappa)































#.............................
# get prev rasters for individual levels
#..............................
# polybb <- osmdata::getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long","lat")

pr$prevrasters <- map(pr$data, 
                      fit_pred_spMLE, outcome = "logitplsmdprev", covar = "1", 
                      long_var = "longnum", lat_var = "latnum",
                      grid.pred = grid.pred, kappa = 0.5, 
                      pred.reps = 10)

pr$prevrasterspred <- purrr::map(pr$prevrasters, "pred")



#.............................
# plot prev rasters
#..............................
prevmaprasterplots <- lapply(pr$prevrasterspred,
                             prevmaprasterplotter, smoothfct = rep(7,3), alpha = 0.5)


#----------------------------------------------------------------------------------------------------
# Save Objects & Write out
#----------------------------------------------------------------------------------------------------

saveRDS(pr, file = "data/derived_data/basic_prevmap_mapping_data.rds")
saveRDS(prevmaprasterplots, file = "results/prevmap_raster_plots.rds")


