# # plot gg
# ggplot() +
#   geom_sf(data = datg$osm_polygons, fill = '#9ecae1', colour = '#9ecae1', alpha = 0.5) +
#   theme_bw()
# 
# jpeg("~/Desktop/test.jpg", width = 8, height = 8, res=200, units = "in")
# PlotObj
# graphics.off()
# 
# 
# 
# 
# datg <- opq (getbb("Democratic Republic of the Congo", featuretype = "country"), memsize = 1e9) %>%
#   add_osm_feature(key = "waterway", value = "river") %>%
#   osmdata_sf (quiet = FALSE)
# 
# 
# bb_poly <- getbb ("Democratic Republic of the Congo", featuretype = "country", format_out = 'polygon')
# 
# datg <- trim_osmdata (datg, bb_poly, exclude = FALSE)
# 
# 




# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
library(sf)
library(ggspatial)

DRCprov <- readRDS("~/Documents/GitHub/VivID_Epi/data/gadm_drclvl1.rds")
colnames(DRCprov) <- tolower(colnames(DRCprov))
colnames(DRCprov)[4] <- "adm1name" # to match the DHS province names
# need to strip accent marks also to match the DHS province names
# https://stackoverflow.com/questions/20495598/replace-accented-characters-in-r-with-non-accented-counterpart-utf-8-encoding
# thanks to @Thomas for this great trick

unwanted_array = list(   'Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                         'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                         'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'ß'='Ss', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                         'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                         'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )

DRCprov$adm1name <- chartr(paste(names(unwanted_array), collapse=''),
                           paste(unwanted_array, collapse=''),
                           DRCprov$adm1name)

bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
dem.raster <- elevatr::get_elev_raster(sf::as_Spatial(DRCprov), z=5) # elevatr expects sp 

# raster to ggplot for color
dem.raster <- raster::crop(dem.raster, sf::as_Spatial(DRCprov), snap='in')
dem.m  <-  raster::rasterToPoints(dem.raster)
dem.df <-  data.frame(lon = dem.m[,1], lat = dem.m[,2], alt = dem.m[,3])

# raster to ggplot for hill share
slope.raster <- raster::terrain(dem.raster, opt='slope')
aspect.raster <- raster::terrain(dem.raster, opt='aspect')
hill.raster <- raster::hillShade(slope.raster, aspect.raster, 40, 270)

hill.m <- raster::rasterToPoints(hill.raster)
hill.df <-  data.frame(lon = hill.m[,1], lat = hill.m[,2], hill = hill.m[,3])


set_overpass_url('https://lz4.overpass-api.de/api/interpreter') # https://github.com/ropensci/osmdata/issues/126
# https://wiki.openstreetmap.org/wiki/Map_Features#Highway
trunkroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "highway", value = "trunk") %>% # The most important roads in a country's system that aren't motorways. (Need not necessarily be a divided highway.) 
  osmdata::osmdata_sf() %>% 
  trim_osmdata(polybb) 
trunkroadsosm <- trunkroadsosm$osm_lines

riverosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "waterway", value = "river") %>% # The linear flow of a river, in flow direction.
  #  add_osm_feature(key = 'name', value = 'Congo', value_exact = FALSE) %>%
 osmdata::osmdata_sf() 
majriver <- riverosm$osm_lines[which(tolower(riverosm$osm_lines$name) %in% c("congo", "ubangi")), ]



plotObj <- ggplot() +
  geom_raster(data=hill.df, aes(lon, lat, fill=hill), alpha=0.5) +
  scale_fill_gradientn(colours = grey.colors(100)) +
  geom_raster(data = dem.df, aes(lon, lat, fill = alt), alpha = .45) +
  scale_fill_gradientn(colours = terrain.colors(100)) +
  geom_sf(data = majriver, colour = '#9ecae1', size = 2, alpha = 0.7) +
  geom_sf(data=DRCprov, fill = NA) +
  theme_bw()

plotObj
