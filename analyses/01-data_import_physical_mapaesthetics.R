#----------------------------------------------------------------------------------------------------
# Purpose of this script is to download map features that will be need for later plotting
# Will use this to mostly make "pretty" maps 
# Am going to recycle the DEM
# http://www.francescobailo.net/2018/08/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(sf)
library(raster)
library(ggspatial)
library(elevatr)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")


#---------------------------------------------------------------------------------
# pull down DRC maps from GADM
#---------------------------------------------------------------------------------
#spatial from GADM -- these are polygon files, doing this is legacy as raster does this nicely...but already fixed naming issue here 
if(!dir.exists(paste0(getwd(), "/data/map_bases/gadm/"))){dir.create(paste0(getwd(), "/data/map_bases/gadm/"), recursive = T)}
DRCprov <- sf::st_as_sf( raster::getData(name = "GADM", country = "CD", level = 1, path = "data/map_bases/gadm/") )
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


#..............................
# Pull down the border cntrs w/ raster
#..............................
brdrcnt <- lapply(c("UGA", "SSD", "CAF", "COG", "COD", "AGO", "ZMB", "TZA", "RWA", "BDI", "GAB", "CMR", "GNQ"), 
                  function(x){
                    ret <- raster::getData(name = "GADM", country = x, level = 0, path = "data/map_bases/gadm/")
                    ret <- sf::st_as_sf(ret)
                    return(ret)
                    
                  })


#..............................
# Pull down ocean
#..............................
# https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip
oceans <- sf::st_read("~/Documents/GitHub/VivID_Epi/data/map_bases/ne_10m_ocean/ne_10m_ocean.shp")

#..............................
# Pull down terrain and hill shading
#..............................
dem.raster <- elevatr::get_elev_raster(sf::as_Spatial(DRCprov), z=5) # elevatr expects sp 

# raster to ggplot for color
dem.raster <- raster::crop(dem.raster, sf::as_Spatial(DRCprov), snap='out')
dem.m  <-  raster::rasterToPoints(dem.raster)
dem.df <-  data.frame(lon = dem.m[,1], lat = dem.m[,2], alt = dem.m[,3])

# raster to ggplot for hill shapes
slope.raster <- raster::terrain(dem.raster, opt='slope')
aspect.raster <- raster::terrain(dem.raster, opt='aspect')
hill.raster <- raster::hillShade(slope.raster, aspect.raster, 40, 270, normalize = T)
hill.m <- raster::rasterToPoints(hill.raster)
hill.df <-  data.frame(lon = hill.m[,1], lat = hill.m[,2], hill = hill.m[,3])


#---------------------------------------------------------------------------------
# write out lists of map bases for later plotting
#---------------------------------------------------------------------------------
prettybasemap_terraincolors <- list(
  geom_raster(data=hill.df, aes(lon, lat, fill=hill)),
  geom_raster(data = dem.df, aes(lon, lat, fill = alt), alpha = 0.7),
  scale_fill_gradientn(colours = terrain.colors(100), guide = F),
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_y = unit(1.25, "cm")),
  vivid_theme,
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)


prettybasemap_nodrc <- list(
  # geom_raster(data=hill.df, aes(lon, lat, fill=hill)) +
  # geom_raster(data = dem.df, aes(lon, lat, fill = alt), alpha = 0.7) +
  #  scale_fill_manual(values = "#bdbdbd", guide = F) +
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  # geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_y = unit(1.25, "cm")),
  vivid_theme,
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)



prettybasemap_nodrc_nonorth <- list(
  # geom_raster(data=hill.df, aes(lon, lat, fill=hill)) +
  # geom_raster(data = dem.df, aes(lon, lat, fill = alt), alpha = 0.7) +
  #  scale_fill_manual(values = "#bdbdbd", guide = F) +
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  # geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  #   ggspatial::annotation_north_arrow(location = "bl", which_north = "true", pad_y = unit(1.25, "cm")),
  vivid_theme,
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)


prettybasemap_hillgrey <- list(
  geom_raster(data=hill.df, aes(lon, lat, fill=hill)),
  # geom_raster(data = dem.df, aes(lon, lat, fill = alt), alpha = 0.7),
  scale_fill_manual(values = "#bdbdbd", guide = F),
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = oceans, fill = "#9ecae1"),
  geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true",
                                    pad_y = unit(1.25, "cm")),
  vivid_theme,
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)




#----------------------------------------------------------------------------------------------------
# Save Objects & Write out
#----------------------------------------------------------------------------------------------------
save(prettybasemap_terraincolors, 
     prettybasemap_hillgrey, 
     prettybasemap_nodrc,
     prettybasemap_nodrc_nonorth,
     file = "data/map_bases/vivid_maps_bases.rda")
saveRDS(DRCprov, file = "data/map_bases/vivid_DRCprov.rds")