#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the night light composite data from VIIRS
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sp)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

#.............................................................................. 
# Night Light Raster Merge
#.............................................................................. 
# night light rasters
N0E60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N060E_v10_c201701311200/SVDNB_npp_20150101-20151231_00N060E_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N0W60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N060W_v10_c201701311200/SVDNB_npp_20150101-20151231_00N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N0W180 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_00N180W_v10_c201701311200/SVDNB_npp_20150101-20151231_00N180W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75E60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N060E_v10_c201701311200/SVDNB_npp_20150101-20151231_75N060E_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75W60 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N060W_v10_c201701311200/SVDNB_npp_20150101-20151231_75N060W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")
N75W180 <- raster::raster("data/raw_data/night_lights/SVDNB_npp_20150101-20151231_75N180W_v10_c201701311200/SVDNB_npp_20150101-20151231_75N180W_vcm-orm-ntl_v10_c201701311200.avg_rade9.tif")

nightlights.rstrs <- list(N0E60, N0W60, N0W180,
                          N75E60, N75W60, N75W180)

nightlights.merge <- Reduce(function(...) merge(...), nightlights.rstrs)
nightlights.drc <- raster::crop(x = nightlights.merge, y = caf)
nightlights.drc <- raster::projectRaster(from = nightlights.drc, to = nightlights.drc,
                                         crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m

# need to write out merge object because it is not saving correctly as .RDS potentially due to layering
# save out this surface
raster::writeRaster(nightlights.drc,  filename = "data/derived_data/vividepi_nightlights_surface.grd ")



# look at data for clusters
summary(values(nightlights.drc))
hist( values(nightlights.drc) )
hist( values(nightlights.drc)[values(nightlights.drc) > 0] )
sum(values(nightlights.drc) < 0) # there are 2237 values less than 0 out of 44920800
# apparently these <0 values can result from too much correction https://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=5888
# given that there are so few, I am going to ignore them

#.............................................................................. 
# Extract Night Light By Cluster
#.............................................................................. 
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
ge.nightlights <- ge %>% 
  dplyr::select(c("hv001", "urban_rura", "geometry")) %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))

# loop through and find means
ge.nightlights$nightlightsmean <- NA

for(i in 1:nrow(ge.nightlights)){
  
  ge.nightlights$nightlightsmean[i] <- 
    raster::extract(x = nightlights.drc, # raster doesn't change 
                    y = sf::as_Spatial(ge.nightlights$geometry[i]),
                    buffer = ge.nightlights$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

# given that the vast majority of these are 0-values
# going to do a zero truncated scale
nonzeroes <- which( ge.nightlights$nightlightsmean > 0 )
ge.nightlights$nightlightsmean_cont_scale_clst <- ge.nightlights$nightlightsmean
ge.nightlights$nightlightsmean_cont_scale_clst[nonzeroes] <- my.scale(ge.nightlights$nightlightsmean[nonzeroes])

# save out
saveRDS(object = ge.nightlights, 
        file = "data/derived_data/vividepi_night_clstmeans.rds")




