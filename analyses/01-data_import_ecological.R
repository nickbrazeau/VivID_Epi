#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(sf)

# set boundaries
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
ape <- sf::read_sf("data/raw_data/redlist_species_data_primate/data_0.shp")
sf::st_crs(ape)
ape <- ape[grepl("pan paniscus|pan troglodytes|gorilla", tolower(ape$BINOMIAL)), ] %>%  # pan trog, pan panisus, gorilla sp
  dplyr::rename(species = BINOMIAL)
drc_ape <- sf::st_crop(x = ape, y = sf::st_as_sf(caf))


#---------------------------------------------------------------------------------
# Precipation Data
#---------------------------------------------------------------------------------
dir.create(path = "data/raw_data/weather_data/CHIRPS/", recursive = T)
heavyRain::getCHIRPS(region = "africa",
                     format = "tifs",
                     tres = "monthly", 
                     sres = 0.05, 
                     begin = as.Date("2013-01-01"),
                     end = as.Date("2014-12-31"),
                     dsn = "data/raw_data/weather_data/CHIRPS/",
                     overwrite = T)

system('gunzip data/raw_data/weather_data/CHIRPS/*')


#---------------------------------------------------------------------------------
# Temperature Data
#---------------------------------------------------------------------------------
# Manually downloaded from LAADS by requesting on their server and using `wget`


##################################################################################
##########                      HOTOSM DATA                       ################
##################################################################################
# read in GE as import
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 
sf::st_crs(ge)


#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp") %>% 
  dplyr::filter(waterway == "river") %>% 
  dplyr::rename(watertype = waterway) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))

wtrply <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp") %>% 
  dplyr::filter(water == "lake") %>% 
  dplyr::rename(watertype = water) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))

wtr <- sf::st_combine(rbind(wtrlns, wtrply))
wtr <-  sf::st_union( wtr )
#......................
# quick look
#......................
sf::st_crs(wtr)
# plot
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
ggplot() + 
  geom_sf(data = sf::st_as_sf(DRCprov)) +
  geom_sf(data = wtr, color = "blue")

#......................
# GC dist
#......................
wtrdist <- sf::st_distance(x = ge,
                           y = wtr,
                           which = "Great Circle")


wtrdist_out <- data.frame(
  hv001 = ge$dhsclust,
  wtrdist_cont_clst = apply(wtrdist, 1, min)
)

# check the results
wtrdist_ge <- wtrdist_out %>% 
  dplyr::rename(dhsclust = hv001) %>% 
  dplyr::left_join(., ge, by = "dhsclust")
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
ggplot() + 
  geom_sf(data = sf::st_as_sf(DRCprov)) +
  geom_sf(data = wtr, color = "blue") +
  geom_point(data = wtrdist_ge, aes(x = longnum, y = latnum, color = wtrdist_cont_clst)) +
  scale_color_viridis_c("Wtr Dist")




#----------------------------------------------------------------------------------------------------
# write out
#----------------------------------------------------------------------------------------------------
dir.create(path = "data/derived_data/", recursive = T)
saveRDS(drc_ape, file = "data/derived_data/drc_ape.rds")
saveRDS(object = wtrdist_out, file = "data/derived_data/hotosm_waterways_dist.rds")
save(wtrply, wtrlns, file = "data/derived_data/hotosm_waterways.RDA")
