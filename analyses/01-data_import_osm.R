#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import GIS variables and perform basic GIS calculations (e.g. various distance metrics)
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)


#---------------------------------------------------------------------------------
# pull down OSM maps
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata
bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')

set_overpass_url('https://lz4.overpass-api.de/api/interpreter') # https://github.com/ropensci/osmdata/issues/126
# https://wiki.openstreetmap.org/wiki/Map_Features#Highway
trunkroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "highway", value = "trunk") %>% # The most important roads in a country's system that aren't motorways. (Need not necessarily be a divided highway.) 
  osmdata::osmdata_sf() %>% 
  trim_osmdata(polybb) 
trunkroadsosm <- trunkroadsosm$osm_lines

primaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "highway", value = "primary") %>% # The next most important roads in a country's system. (Often link larger towns.)
  osmdata::osmdata_sf() %>% 
  trim_osmdata(polybb) 
primaryroadsosm <- primaryroadsosm$osm_lines

secondaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
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
# 
# hospitalosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "amenity", value = "hospital") %>% 
#   osmdata::osmdata_sf() %>% 
#   trim_osmdata(polybb) 
# hospitalosm <- hospitalosm$osm_points
# 
# docosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "amenity", value = "doctor") %>% # A doctor's practice / surgery.
#   osmdata::osmdata_sf() %>% 
#   trim_osmdata(polybb) 
# docosm <- docosm$osm_points

riverosm <- osmdata::opq(bbox = bb, memsize = 1e11 ) %>%
  add_osm_feature(key = "waterway", value = "river") %>% # The linear flow of a river, in flow direction.
#  add_osm_feature(key = 'name', value = 'Congo', value_exact = FALSE) %>%
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb, exclude = F)

majriver <- riverosm$osm_lines[which(tolower(riverosm$osm_lines$name) %in% c("congo", "ubangi")), ]



# riverbankosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "waterway", value = "riverbank") %>% # A wide river as defined by its area.
#   trim_osmdata(polybb) %>% 
#   osmdata::osmdata_sf()

# streamosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "waterway", value = "stream") %>% # The linear flow of a river, in flow direction.
#   trim_osmdata(polybb) %>% 
#   osmdata::osmdata_sf()

# lakeosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "water", value = "lake") %>% # Any body of water, from natural such as a lake or pond to artificial like moat or canal
#   osmdata::osmdata_sf() %>% 
#   trim_osmdata(polybb, exclude = F) 
# lakeosm <- lakeosm$osm_multipolygons



