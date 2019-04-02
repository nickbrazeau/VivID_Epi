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

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

mp$data <- lapply(mp$data, function(x){
  return( dplyr::left_join(x = x, y = ge, by ="hv001") )
})


#......................
# Summarize by Cluster
#......................



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



#..........................
# spatial kernel densities
#..........................
# cases <- dt[dt$pv18s == 1, ]
# 
# casedens <- MASS::kde2d(x = cases$longnum, 
#                  y = cases$latnum,
#                  n=1e3,
#             lims = c(bb[1,], bb[2,]))
# 
# contour(casedens)
# 
# 
# controls <- dt[dt$pv18s == 0, ]
# condens <- MASS::kde2d(x = controls$longnum, 
#                     y = controls$latnum,
#                     n=1e3,
#                     lims = c(bb[1,], bb[2,]))
# 
# contour(condens)



