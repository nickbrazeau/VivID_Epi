#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# landcoverage data for crops for farmer variable
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("R/00-functions_basic.R")

# create bounding box of Central Africa for mrm
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"

# create mask 
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

#.............................................................................. 
# Read in Land Coverage
#.............................................................................. 
landcov2013 <- raster::raster("data/raw_data/land_coverage/dataset-satellite-land-cover-902f410c-4e52-4fae-badd-2fc35ec86c62/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2013-v2.0.7cds.nc")
# sanity
raster::compareCRS(landcov2013, raster::crs("+init=epsg:4326"))

# speed
landcov2013.drc <- raster::crop(x = landcov2013, y = caf)
# mask out non DRC 
landcov2013.drc <- raster::mask(landcov2013.drc, DRC)

#.............................................................................. 
# Lift Over to Binary Cropland yes or no
#.............................................................................. 
# using data dictionary from here: https://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-QuickUserGuide-LC-Maps_v2-0-7.pdf
origvals <- raster::values(landcov2013.drc)
newvals <- ifelse(origvals %in% c(10, 20, 30, 40), 1, 
                  ifelse(is.na(origvals), NA, 0))
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
ge.croopland$cropprop_cont_scale_clst <- my.scale(logit(ge.croopland$cropprop, tol = 1e-3))

# drop extraneous geometry
sf::st_geometry(ge.croopland) <- NULL


# save out
saveRDS(object = ge.croopland, 
        file = "data/derived_data/vividepi_cropland_propmeans.rds")





