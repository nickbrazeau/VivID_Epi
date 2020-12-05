#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle distance to public Health Sites
# in the DRC as a proxy to healthcare accessibility
# recent work has made this available at an even higher resolution: 
#  https://malariaatlas.org/research-project/accessibility_to_healthcare/
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(sp)
library(sf)

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"


#..................................
# import data
#..................................
# Get cluster locations
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# drop observations with missing geospatial data 
ge <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::select(c("hv001", "longnum", "latnum", "urban_rura")) %>% 
  dplyr::mutate(urban_rura = as.character(urban_rura)) %>% # coerce to char so attr list doesn't mess with dplyr
  dplyr::filter(!duplicated(.))
# sanity check
sf::st_crs(ge)
identicalCRS(sf::as_Spatial(ge), caf)
# liftover to conform with rgdal updates http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
ge <- sp::spTransform(sf::as_Spatial(ge), CRSobj = sp::CRS("+init=epsg:4326"))
identicalCRS(ge, caf)
# back to tidy 
ge <- sf::st_as_sf(ge)



#............................................................
# get health care calculated distances
#...........................................................
hlthdist_walk <- raster::raster("data/raw_data/hlthdist/2020_walking_only_travel_time_to_healthcare.geotiff")
hlthdist_road <- raster::raster("data/raw_data/hlthdist/2020_motorized_travel_time_to_healthcare.geotiff")

hlthdist_stack <- raster::stack(hlthdist_walk, hlthdist_road)
hlthdist_stack <- raster::crop(x = hlthdist_stack, y = caf)
hlthdist <- raster::calc(hlthdist_stack, mean, na.rm = T)

# sanity check
sf::st_crs(hlthdist)
# crop for speed

# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
hlthdist <- raster::mask(x = hlthdist, mask = DRCprov)


#............................................................
# new get buffer around urban versus rural cluster
#...........................................................
# hlthdist.mean <- ge[,c("hv001", "geometry", "urban_rura")] %>% 
#   dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))
# hlthdist.mean <- hlthdist.mean[!duplicated(hlthdist.mean$hv001),]

hlthdist.mean <- ge[,c("hv001", "geometry", "urban_rura")] 
hlthdist.mean <- hlthdist.mean[!duplicated(hlthdist.mean$hv001),]



# Drop in a for loop again to account for DHS buffering
# note the 0.05 degree resolution is approximately 6km, so the buffer
# for urbanicity shouldn't be doing anything... but to be consistent with 
# DHS "The Geospatial Covariate Datasets Manual", we will do it
hlthdist.mean$hlthdist_cont_clst <- NA

for(i in 1:nrow(hlthdist.mean)){
  hlthdist.mean$hlthdist_cont_clst[i] <- 
    raster::extract(x = hlthdist, 
                    y = sf::as_Spatial(hlthdist.mean$geometry[i]),
                    buffer = 100e3, # 100 km regardless of urban/rural
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}


#......................
# inspect results
#......................
# aggregate quickly so we don't overwhelm ggplot
hlthdist.agg <- raster::aggregate(hlthdist, fact = 10, fun = mean)

hlthdist.mean %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2]) %>% 
  ggplot() + 
  geom_sf(data = sf::st_as_sf(DRCprov)) +
  ggspatial::layer_spatial(data = hlthdist.agg,
                           aes(fill = stat(band1)),
                           alpha = 0.9,
                           na.rm = T) +
  geom_point(aes(x = longnum, y = latnum, color = hlthdist_cont_clst)) +
  scale_fill_viridis_c(option="plasma", direction = 1) +
  scale_color_viridis_c(option="viridis")



#..........
# Now write out
#..........
dir.create("data/derived_data/", recursive = TRUE)
saveRDS(object = hlthdist.mean, file = "data/derived_data/hlthdist_out_wlk_trvltime.rds")
