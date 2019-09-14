library(tidyverse)
library(raster)
library(sp)
library(sf)
#----------------------------------------------------------------------------------------------------
# Read Raster Temp and Precip
#----------------------------------------------------------------------------------------------------
readRasterBB <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = bb)
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m
  
  return(ret)
  
}

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"


#.............................................................
# Temperature (MODIS/LAADS) 
#.............................................................
# FILENAME DECONVOLUTION: An example of an existing file is below. 
# https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MOD02QKM/2007/018/MOD02QKM.A2007018.0105.006.2014227230926.hdf
# This path should return a MODIS Terra quarter kilometer (250 m) top of atmosphere reflectance product for year 2007, day-of-year 018 (i.e. January 18), from collection 6.
# Note to change values back to celsius: https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius

# for temperature, need to mask water sources which have too low temp readings
load("data/derived_data/hotosm_waterways.RDA")
oceans <- sf::st_read("data/map_bases/ne_10m_ocean/ne_10m_ocean.shp") %>% 
  sf::st_crop(caf)

tempfiles <- list.files(path = "data/raw_data/weather_data/LAADS_NASA/", full.names = T,
                        pattern = "_Night_CMG.tif")
tempfrst <- lapply(tempfiles, readRasterBB, bb = caf)


# rescale values to celsius
tempfrst <- lapply(tempfrst, function(x){
  # mask water
  x <- raster::mask(x, sf::as_Spatial(wtrlns))
  x <- raster::mask(x, sf::as_Spatial(wtrply))
  x <- raster::mask(x, sf::as_Spatial(oceans))
  
  
  vals <- raster::values(x) 
  vals <- ifelse(vals <= 7500, NA, vals) # improper values
  vals <- (vals * 0.02) - 273.15
  raster::values(x) <- vals
  return(x)
})

saveRDS(tempfrst, file = "data/raw_data/weather_data/LAADS_NASA/MODIS_temp_masked.RDS")
saveRDS(tempfrst, file = "MODIS_temp_masked.RDS")
cat("Finished successfully")

