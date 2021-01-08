#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the night light composite data from VIIRS
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(malariaAtlas)
library(raster)
library(sp)
library(sf)
source("R/00-functions_basic.R")


# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
sf::st_crs(DRCprov)
#.............................................................................................. 
##### Access to City ######
#.............................................................................................. 
#.............................................................................. 
# access surface of accessibility to cities from malaria atlas project
#.............................................................................. 
# find access raster
available_rasters <- malariaAtlas::listRaster()
dir.create("data/raw_data/MAPrasters/", recursive = TRUE)
# create bounding box of Central Africa for download speed
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"

# note, file is time stamped on my computer (would be different from a new user)
if(!file.exists("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_11_22.tiff")) {
  # download from MAP
  malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                          shp = caf,
                          file_path = "data/raw_data/MAPrasters/") 
}


# read in 
accraw <- raster::raster("data/raw_data/MAPrasters/getRaster/2015_accessibility_to_cities_v1.0_latest_10_.18_40_8_2020_11_22.tiff")
raster::crs(accraw)
identicalCRS(accraw, caf)
identicalCRS(accraw, sf::as_Spatial(DRCprov))

# crop for speed
accraw <- raster::crop(x = accraw, y = caf)
# mask
access.drc <- raster::mask(accraw, DRCprov)

# tidy up 
summary(values(access.drc))
head(table(values(access.drc)))
values(access.drc)[values(access.drc) < 0] <- NA
summary(values(access.drc))
hist( values(access.drc) )

# save out this surface
saveRDS(object = access.drc, file = "data/derived_data/vividepi_accessurban_surface.rds")



#.............................................................................. 
# Extract access By Cluster
#.............................................................................. 
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
ge.accurban <- ge %>% 
  dplyr::select(c("hv001", "urban_rura", "geometry")) %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))

# loop through and find means
ge.accurban$accmean <- NA

for(i in 1:nrow(ge.accurban)){
  
  ge.accurban$accmean[i] <- 
    raster::extract(x = access.drc, # raster doesn't change 
                    y = sf::as_Spatial(ge.accurban$geometry[i]),
                    buffer = ge.accurban$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

# distribution looks approximately logged, which is expected 
hist(ge.accurban$accmean)
ge.accurban$accmean_cont_scale_clst <- my.scale(log(ge.accurban$accmean + 0.5)) # offset for zeroes
hist(ge.accurban$accmean_cont_scale_clst)

# drop extraneous geometry
sf::st_geometry(ge.accurban) <- NULL

# save out
saveRDS(object = ge.accurban, 
        file = "data/derived_data/vividepi_accurban_clstmeans.rds")




#.............................................................................................. 
##### Health Distance ######
#.............................................................................................. 
#............................................................
# get health care calculated distances
# just looking at walking distance
#...........................................................
hlthdist <- raster::raster("data/raw_data/hlthdist/2020_walking_only_travel_time_to_healthcare.geotiff")

# sanity check
sf::st_crs(hlthdist)
identicalCRS(hlthdist, caf)
identicalCRS(hlthdist, sf::as_Spatial(DRCprov))

# crop for speed
hlthdist <- raster::crop(x = hlthdist, y = caf)
# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
hlthdist.drc <- raster::mask(x = hlthdist, mask = DRCprov)



#.............................................................................. 
# Extract health distance By Cluster
#.............................................................................. 
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
ge.hlthdisturban <- ge %>% 
  dplyr::select(c("hv001", "urban_rura", "geometry")) %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))

# loop through and find means
ge.hlthdisturban$hlthdist <- NA

for(i in 1:nrow(ge.hlthdisturban)){
  
  ge.hlthdisturban$hlthdist[i] <- 
    raster::extract(x = hlthdist.drc, # raster doesn't change 
                    y = sf::as_Spatial(ge.hlthdisturban$geometry[i]),
                    buffer = ge.hlthdisturban$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

# distribution looks approximately logged, which is expected since it is a distance time
summary(ge.hlthdisturban$hlthdist)
hist(ge.hlthdisturban$hlthdist)
ge.hlthdisturban$hlthdist_cont_scale_clst <- my.scale(log(ge.hlthdisturban$hlthdist + 0.5)) # offset for zeroes
hist(ge.hlthdisturban$hlthdist_cont_scale_clst)

# drop extraneous geometry
sf::st_geometry(ge.hlthdisturban) <- NULL

# save out
saveRDS(object = ge.hlthdisturban, 
        file = "data/derived_data/vividepi_hlthdist_clstmeans.rds")


#............................................................
# Access to Citie and Heaalth Distance Correlations
#...........................................................
cor(values(access.drc),
    values(hlthdist.drc),
    use = "na.or.complete")
