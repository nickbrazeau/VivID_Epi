#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the night light composite data from MCD12Q1 2013
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3

# create bounding box of Central Africa for mrm
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

#.............................................................................. 
# Read in Land Coverage
#.............................................................................. 
landcov2013 <- raster::raster("data/raw_data/land_coverage/dataset-satellite-land-cover-902f410c-4e52-4fae-badd-2fc35ec86c62/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2013-v2.0.7cds.nc")
landcov2013.drc <- raster::crop(x = landcov2013, y = caf)
# mask out non DRC 
landcov2013.drc <- raster::mask(landcov2013.drc, DRCprov)
landcov2013.drc <- raster::projectRaster(from = landcov2013.drc, to = landcov2013.drc,
                                         crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m

#.............................................................................. 
# Lift Over to Binary Cropland yes or no
#.............................................................................. 
# using data dictionary from here: https://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-QuickUserGuide-LC-Maps_v2-0-7.pdf
origvals <- raster::values(landcov2013.drc)
newvals <- ifelse(origvals %in% c(10, 20, 30, 40), 1, 0)
raster::values(landcov2013.drc) <- newvals

saveRDS(object = landcov2013.drc, file = "data/derived_data/vividepi_cropland_surface.rds")

#.............................................................................. 
# Extract Crop Land Cluster
#.............................................................................. 
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
ge.croopland <- ge %>% 
  dplyr::select(c("hv001", "urban_rura", "geometry")) %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))

# loop through and find means
ge.croopland$cropprop <- NA

for(i in 1:nrow(ge.croopland)){
  
  ge.croopland$cropprop[i] <- 
    raster::extract(x = landcov2013.drc, # this doesn't change
                    y = sf::as_Spatial(ge.croopland$geometry[i]),
                    buffer = ge.croopland$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

# liftover
ge.croopland$cropprop_cont_scale_clst <- my.scale(logit(ge.croopland$cropprop, tol = tol))

# save out
saveRDS(object = ge.croopland, 
        file = "data/derived_data/vividepi_cropland_propmeans.rds")





