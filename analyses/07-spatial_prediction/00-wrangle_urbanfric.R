#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the night light composite data from VIIRS
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
remotes::install_github("malaria-atlas-project/malariaAtlas")
library(malariaAtlas)
library(raster)
library(sp)
library(sf)
source("R/00-functions_basic.R")
tol <- 1e-3


# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
sf::st_crs(DRCprov)
#.............................................................................. 
# Friction Surface from malaria atlas project
#.............................................................................. 
# find friction raster
available_rasters <- malariaAtlas::listRaster()
dir.create("data/raw_data/MAPrasters/", recursive = TRUE)
# create bounding box of Central Africa for download speed
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"

# download from MAP
malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                        shp = caf,
                        file_path = "data/raw_data/MAPrasters/") 

# read in 
fricraw <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_11_22.tiff")
raster::crs(fricraw)
# sanity
fricraw <- raster::projectRaster(from = fricraw, to = fricraw,
                                 crs = sf::st_crs("+init=epsg:4326")) 
friction.drc <- raster::mask(fricraw, DRCprov)

# tidy up 
summary(values(friction.drc))
table(values(friction.drc))
values(friction.drc)[values(friction.drc) < 0] <- NA
summary(values(friction.drc))
hist( values(friction.drc) )

# save out this surface
saveRDS(object = friction.drc, file = "data/derived_data/vividepi_frictionurban_surface.rds")



#.............................................................................. 
# Extract Friction By Cluster
#.............................................................................. 
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
ge.fricurban <- ge %>% 
  dplyr::select(c("hv001", "urban_rura", "geometry")) %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))

# loop through and find means
ge.fricurban$fricmean <- NA

for(i in 1:nrow(ge.fricurban)){
  
  ge.fricurban$fricmean[i] <- 
    raster::extract(x = friction.drc, # raster doesn't change 
                    y = sf::as_Spatial(ge.fricurban$geometry[i]),
                    buffer = ge.fricurban$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

# distribution looks approximately logged, which is expected 
hist(ge.fricurban$fricmean)
ge.fricurban$frctmean_cont_scale_clst <- my.scale(log(ge.fricurban$fricmean + tol))
hist(ge.fricurban$frctmean_cont_scale_clst)

# save out
saveRDS(object = ge.fricurban, 
        file = "data/derived_data/vividepi_fricurban_clstmeans.rds")




