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
load("data/map_bases/vivid_maps_bases.rda")

# Summarize by Cluster
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
mp <- mp %>% 
  dplyr::filter(maplvl == "hv001") 


#----------------------------------------------------------------------------------------------------
# Smoothed Guassian Maps
#----------------------------------------------------------------------------------------------------

# bind those to a tibble
pr <- dplyr::bind_rows(mp$data) %>% 
  dplyr::group_by(plsmdmspec) %>% 
  tidyr::nest()


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
prevmaprasterplots <- map(prevmaprasterplots, function(x){return(x + prettybasemap_nodrc)})




#----------------------------------------------------------------------------------------------------
# Save Objects & Write out
#----------------------------------------------------------------------------------------------------

saveRDS(pr, file = "data/derived_data/basic_prevmap_mapping_data.rds")
saveRDS(prevmaprasterplots, file = "results/prevmap_raster_plots.rds")


