#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the covars needed to fit the prov carbayes models
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sf)
source("R/00-functions_basic.R")
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")

#-------------------------------------------------------------------------
# Aggregate Covariates
# means within each adm1
#-------------------------------------------------------------------------
precip <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_precip_study_period_effsurface.rds")
temp <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_temperature_study_period_effsurface.rds")
cropland <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_cropland_surface.rds")
nightlights <- raster::raster("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_nightlights_surface.grd")


extract_agg_raster_polygon <- function(rstrlyr, plygn){
  vals <- raster::extract(x = rstrlyr, y = sf::as_Spatial(plygn),
                          fun = mean,
                          na.rm = T, 
                          sp = F
  )
  return(as.vector(vals))
  
}

# split admins
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
adm1 <- mp$data[mp$plsmdmspec == "pv18s" & mp$maplvl == "adm1name"][[1]] %>% # doesn't matter which prev we use
  dplyr::select(c("adm1name", "geometry")) %>% 
  dplyr::filter(!duplicated(.))

pvcovar <- adm1 %>% 
  dplyr::select(c("adm1name"))

# split into list
adm1 <- split(adm1, 1:nrow(adm1))



pvcovar$precip <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = precip) )
pvcovar$temp <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = temp) )
pvcovar$crop <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = cropland) )
pvcovar$nightlight <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = nightlights) )


pvcovar$precip_scale <- my.scale(pvcovar$precip)
pvcovar$temp_scale <- my.scale(pvcovar$temp)
pvcovar$crop_scale <- my.scale(logit(pvcovar$crop, tol = 1e-3))
pvcovar$nightlight_scale <- my.scale(pvcovar$nightlight)
sf::st_geometry(pvcovar) <- NULL


#..............................................................
# save out
#..............................................................
saveRDS(pvcovar, file = "data/derived_data/vividepi_prov_covars_bayesian_fit.RDS")
