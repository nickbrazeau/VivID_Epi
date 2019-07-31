#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(heavyRain)
library(sf)



#---------------------------------------------------------------------------------
# Precipation Data
#---------------------------------------------------------------------------------
dir.create("data/raw_data/weather_data/CHIRPS/", recursive = T)
heavyRain::getCHIRPS(region = "africa",
                     format = "tifs",
                     tres = "monthly", 
                     sres = 0.05, # near same resolution as Manny pulled down
                     begin = as.Date("2013-01-01"),
                     end = as.Date("2014-12-31"),
                     dsn = "data/raw_data/weather_data/CHIRPS/",
                     overwrite = T)

system('gunzip data/raw_data/weather_data/CHIRPS/*')




#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
ape <- sf::read_sf("data/raw_data/redlist_species_data_primate/data_0.shp")
ape <- ape[grepl("pan paniscus|pan troglodytes|gorilla", tolower(ape$BINOMIAL)), ] %>%  # pan trog, pan panisus, gorilla sp
  dplyr::rename(species = BINOMIAL)


bb <- osmdata::getbb("Democratic Republic of the Congo", 
                     featuretype = "country",
                     format_out = "sf_polygon")
drc_ape <- sf::st_crop(x = ape, y = bb)


saveRDS(drc_ape, file = "data/redlist_species_data_primate/drc_ape.rds")






##################################################################################
##########                      HOTOSM DATA                       ################
##################################################################################
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 

# Note manually downloading these from site
# https://data.humdata.org/dataset/

#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp") %>% 
  dplyr::filter(waterway %in% c("stream", "river", "riverbank")) %>% 
  dplyr::rename(water = waterway) %>% 
  dplyr::select(c("osm_id", "water", "geometry"))

wtrply <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp") %>% 
  dplyr::filter(water == "lake") %>% 
  dplyr::select(c("osm_id", "water", "geometry"))

wtr <- sf::st_combine(rbind(wtrlns, wtrply))
wtr <-  sf::st_union( wtr )

wtrdist <- sf::st_distance(x = ge,
                           y = wtr)
wtrdist_out <- data.frame(
  hv001 = ge$dhsclust,
  wtrdist_cont_clst = apply(wtrdist, 1, min)
)


#----------------------------------------------------------------------------------------------------
# Health Sites
#----------------------------------------------------------------------------------------------------
hlthsites <- sf::read_sf("data/raw_data/hotosm_data/hotosm_drc_healthsites_shapefiles/healthsites.shp") 
htlhdist <- sf::st_distance(x = ge,
                            y = hlthsites)

hlthdist_out <- data.frame(
  hv001 = ge$dhsclust,
  hlthdist_cont_clst = apply(htlhdist, 1, min)
)



#----------------------------------------------------------------------------------------------------
# write out
#----------------------------------------------------------------------------------------------------
saveRDS(object = wtrdist_out, file = "data/derived_data/hotosm_waterways_dist.rds")
saveRDS(object = hlthdist_out, file = "data/derived_data/hotosm_healthsites_dist.rds")
