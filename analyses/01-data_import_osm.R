#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data from open street map 
# primarily interested in road and water networks
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
library(rgeos)
library(rgdal)
library(sf)
tol <- 1e-3
#---------------------------------------------------------------------------------
# pull down OSM maps
#---------------------------------------------------------------------------------
api_list <- c('http://overpass-api.de/api/interpreter',
              'https://lz4.overpass-api.de/api/interpreter',
              'https://z.overpass-api.de/api/interpreter',
              'https://overpass.kumi.systems/api/interpreter')




# https://github.com/ropensci/osmdata
bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')

# https://wiki.openstreetmap.org/wiki/Map_Features

#---------------------------------------------------------------------------------
# Road Network
#---------------------------------------------------------------------------------
api_to_use <- sample(1:length(api_list), 1)
set_overpass_url(api_list[api_to_use]) 

trunkroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "highway", value = "trunk") %>% # The most important roads in a country's system that aren't motorways. (Need not necessarily be a divided highway.) 
  osmdata::osmdata_sf() %>% 
  trim_osmdata(polybb) 
trunkroadsosm <- trunkroadsosm$osm_lines

primaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e11 ) %>%
  add_osm_feature(key = "highway", value = "primary") %>% # The next most important roads in a country's system. (Often link larger towns.)
  osmdata::osmdata_sf() %>% 
  trim_osmdata(polybb) 
primaryroadsosm <- primaryroadsosm$osm_lines

secondaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e11 ) %>%
  add_osm_feature(key = "highway", value = "secondary") %>% # The next most important roads in a country's system. (Often link towns.)
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb)
secondaryroadsosm <- secondaryroadsosm$osm_lines

# 
# tertiaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "highway", value = "tertiary") %>% # The next most important roads in a country's system. (Often link smaller towns and villages)
#   osmdata::osmdata_sf() %>% 
#   trim_osmdata(polybb) 
# tertiaryroadsosm <- tertiaryroadsosm$osm_lines

roadnet <- sf::st_union(trunkroadsosm, primaryroadsosm) %>% 
  sf::st_union(., secondaryroadsosm) %>% 
  rgeos::gSimplify(., tol = tol)

#---------------------------------------------------------------------------------
# Hospitals & Clinics & Practices
#---------------------------------------------------------------------------------
api_to_use <- sample(1:length(api_list), 1)
set_overpass_url(api_list[api_to_use]) 

hospitalosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "amenity", value = "hospital") %>%
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb)
hospitalosm <- hospitalosm$osm_points %>% 
  dplyr::select(c("osm_id", "geometry")) %>% 
  dplyr::mutate(level = "hosp")

clinosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "amenity", value = "clinic") %>% # A medium-sized medical facility or health centre.
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb)
clinosm <- clinosm$osm_points %>% 
  dplyr::select(c("osm_id", "geometry")) %>% 
  dplyr::mutate(level = "clin")


docosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "amenity", value = "doctors") %>% # A doctor's practice / surgery.
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb)
docosm <- docosm$osm_points %>% 
  dplyr::select(c("osm_id", "geometry")) %>% 
  dplyr::mutate(level = "dctr")

hlthoffc <- sf::st_union(hospitalosm, clinosm) %>% 
  sf::st_union(., docosm)

#---------------------------------------------------------------------------------
# Bodies of Water 
#---------------------------------------------------------------------------------
api_to_use <- sample(1:length(api_list), 1)
set_overpass_url(api_list[api_to_use]) 

riverosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "water", value = "river") %>% 
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb, exclude = F)
riverout <- sf::st_union(x = riverosm$osm_lines, y = riverosm$osm_polygons) %>% 
  sf::st_union(x = ., y = riverosm$osm_multilines) %>% 
  sf::st_union(x = ., y = riverosm$osm_multipolygons)


streamosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "water", value = "stream") %>% # The linear flow of a river, in flow direction.
  trim_osmdata(polybb) %>%
  osmdata::osmdata_sf()
streamout <- sf::st_union(x= streamosm$osm_lines, y = streamosm$osm_polygons) %>% 
  sf::st_union(x = ., y = streamosm$osm_multilines) %>% 
  sf::st_union(x = ., y = streamosm$osm_multipolygons)

  
lakeosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "water", value = "lake") %>% # Any body of water, from natural such as a lake or pond to artificial like moat or canal
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb, exclude = F)
lakeout <- sf::st_union(x= lakeosm$osm_lines, y = lakeosm$osm_polygons) %>% 
  sf::st_union(x = ., y = lakeosm$osm_multilines) %>% 
  sf::st_union(x = ., y = lakeosm$osm_multipolygons)



wtr <- sf::st_union(riverosm, streamosm) %>% 
  sf::st_union(., lakeosm) %>% 
  rgeos::gSimplify(., tol = tol)


saveRDS(roadnet, file = "data/derived_data/osm_roads.rds")
saveRDS(hlthoffc, file = "data/derived_data/osm_healthoffices.rds")
saveRDS(majrivers, file = "data/derived_data/osm_majrivers.rds")
saveRDS(wtr, file = "data/derived_data/osm_water.rds")


