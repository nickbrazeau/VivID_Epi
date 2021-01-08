#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs wrangling script
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(sf)
library(sp)

# set boundaries
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
# http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
sp::proj4string(caf) <- "+init=epsg:4326"

#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
ape <- sf::read_sf("data/raw_data/redlist_species_data_primate/data_0.shp")
# for sanity
sf::st_crs(ape)
identicalCRS(caf, sf::as_Spatial(ape))
# subset to greater apes 
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
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
# sanity check
sf::st_crs(ge)
identicalCRS(sf::as_Spatial(ge), caf)
# liftover to conform with rgdal updates http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
ge <- sp::spTransform(sf::as_Spatial(ge), CRSobj = sp::CRS("+init=epsg:4326"))
identicalCRS(ge, caf)
# back to tidy 
ge <- sf::st_as_sf(ge)

#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp") %>% 
  dplyr::filter(waterway == "river") %>% 
  dplyr::rename(watertype = waterway) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))
# sanity check
sf::st_crs(wtrlns)

wtrply <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp") %>% 
  dplyr::filter(water == "lake") %>% 
  dplyr::rename(watertype = water) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))
# sanity check
sf::st_crs(wtrply)

#......................
# quick look
#......................
# plot
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
ggplot() + 
  geom_sf(data = sf::st_as_sf(DRCprov)) +
  geom_sf(data = wtrlns, color = "blue") + 
  geom_sf(data = wtrply, color = "blue") + 
  geom_sf(data = ge, color = "red")

#......................
# GC dist
#......................
riverdist <- sf::st_distance(x = ge,
                             y = wtrlns,
                             which = "Great Circle")

lakedist <- sf::st_distance(x = ge,
                            y = wtrply,
                            which = "Great Circle")


wtrdist_out <- data.frame(
  hv001 = ge$hv001,
  river = apply(riverdist, 1, min),
  lake = apply(lakedist, 1, min)) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::mutate(wtrdist_cont_clst = min(river, lake)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c("hv001", "wtrdist_cont_clst"))

# check the results
wtrdist_ge <- wtrdist_out %>% 
  dplyr::left_join(., ge, by = "hv001") %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2])
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
ggplot() + 
  geom_sf(data = sf::st_as_sf(DRCprov)) +
  geom_sf(data = wtrlns, color = "blue") + 
  geom_sf(data = wtrply, color = "blue") + 
  geom_point(data = wtrdist_ge, aes(x = longnum, y = latnum, color = wtrdist_cont_clst)) +
  scale_color_viridis_c("Wtr Dist")




#----------------------------------------------------------------------------------------------------
# write out
#----------------------------------------------------------------------------------------------------
dir.create(path = "data/derived_data/", recursive = T)
saveRDS(drc_ape, file = "data/derived_data/drc_ape.rds")
saveRDS(object = wtrdist_out, file = "data/derived_data/hotosm_waterways_dist.rds")
save(wtrply, wtrlns, file = "data/derived_data/hotosm_waterways.RDA")

