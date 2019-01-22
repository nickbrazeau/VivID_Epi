#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
library(sf)
source("analyses/00-functions_basic.R")
#---------------------------------------------------------------------------------
# pull down DRC maps from GADM
#---------------------------------------------------------------------------------
#spatial from GADM -- these are polygon files, doing this is legacy as raster does this nicely...but already fixed naming issue here 
if(!dir.exists(paste0(getwd(), "/data/gadm/"))){dir.create(paste0(getwd(), "/data/gadm/"), recursive = T)}
drclvl0 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_0_sf.rds", httr::write_disk(path="data/gadm/gadm_drclvl0.rds", overwrite = T))
drclvl1 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_1_sf.rds", httr::write_disk(path="data/gadm/gadm_drclvl1.rds", overwrite = T))
drclvl2 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_2_sf.rds", httr::write_disk(path="data/gadm/gadm_drclvl2.rds", overwrite = T))
drclvl3 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_3_sf.rds", httr::write_disk(path="data/gadm/gadm_drclvl3.rds", overwrite = T))

DRCprov <- readRDS("data/gadm_drclvl1.rds")
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
brdrcnt <- lapply(c("UGA", "SSD", "CAF", "COG", "AGO", "ZMB", "TZA", "RWA", "BDI", "GAB", "CMR", "GNQ"), 
                  function(x){
                    ret <- raster::getData(name = "GADM", country = x, level = 0, path = "data/gadm/")
                    ret <- sf::st_as_sf(ret)
                    return(ret)
                    
                  })

# brdrcnt_comb <- sf::st_union(brdrcnt[[1]], brdrcnt[[2]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[3]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[4]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[5]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[6]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[7]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[8]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[9]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[10]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[11]], by_feature = T) %>% 
#   sf::st_union(., brdrcnt[[12]])
# unfortunately this knocks out borders (logically)... will have to keep as seperate layers
#..............................
# Pull down ocean
#..............................
# https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_ocean.zip
oceans <- sf::st_read("data/ne_10m_ocean/ne_10m_ocean.shp")
#---------------------------------------------------------------------------------
# pull down OSM maps
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata
# http://www.francescobailo.net/2018/08/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/


bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
dem.raster <- elevatr::get_elev_raster(sf::as_Spatial(DRCprov), z=5) # elevatr expects sp 

# raster to ggplot for color
dem.raster <- raster::crop(dem.raster, sf::as_Spatial(DRCprov), snap='out')
dem.m  <-  raster::rasterToPoints(dem.raster)
dem.df <-  data.frame(lon = dem.m[,1], lat = dem.m[,2], alt = dem.m[,3])

# raster to ggplot for hill share
slope.raster <- raster::terrain(dem.raster, opt='slope')
aspect.raster <- raster::terrain(dem.raster, opt='aspect')
hill.raster <- raster::hillShade(slope.raster, aspect.raster, 40, 270, normalize = T)

hill.m <- raster::rasterToPoints(hill.raster)
hill.df <-  data.frame(lon = hill.m[,1], lat = hill.m[,2], hill = hill.m[,3])



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

# secondaryroadsosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
#   add_osm_feature(key = "highway", value = "primary") %>% # The next most important roads in a country's system. (Often link towns.)
#   osmdata::osmdata_sf() %>% 
#   trim_osmdata(polybb) 
# secondaryroadsosm <- secondaryroadsosm$osm_lines
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

riverosm <- osmdata::opq(bbox = bb, memsize = 1e9 ) %>%
  add_osm_feature(key = "waterway", value = "river") %>% # The linear flow of a river, in flow direction.
#  add_osm_feature(key = 'name', value = 'Congo', value_exact = FALSE) %>%
  osmdata::osmdata_sf() %>%
  trim_osmdata(polybb, exclude = F)

majriver <- riverosm$osm_lines[which(tolower(riverosm$osm_lines$name) %in% c("congo", "ubangi")), ]
# 


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




#---------------------------------------------------------------------------------
# Pull Down Stamen Maps (esp terrain)
#---------------------------------------------------------------------------------

drc_stamen_back_terrain <- ggmap::get_stamenmap(bb, zoom = 5, 
                                                maptype = "terrain-background")

#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
primate <- sf::read_sf("data/redlist_species_data_primate/data_0.shp")
ape <- primate[grepl("pan paniscus|pan troglodytes|gorilla", tolower(primate$BINOMIAL)), ] # pan trog, pan panisus, gorilla sp


#---------------------------------------------------------------------------------
# Climate Data from Outside Sources
#---------------------------------------------------------------------------------

# From OJ (via slack on 1/10/2019)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056487
# https://www.nature.com/articles/ncomms1879#supplementary-information
# https://betterexplained.com/articles/an-interactive-guide-to-the-fourier-transform/
load("data/imperial_share/admin_units_seasonal.rda")
admin_units_seasonal <- admin_units_seasonal %>% 
  filter(country == "Democratic Republic of the Congo")
# this is from 2010 and just the 13 DRC prov... 

# if this doesn't pan out, check the GSODR package from Ropensci -- https://cran.r-project.org/web/packages/GSODR/index.html


#---------------------------------------------------------------------------------
# Climate Data from DHS, this is a rich dataset ... will need to think/discuss as group
#---------------------------------------------------------------------------------
gc <- readr::read_csv("datasets/CDGC62FL/CDGC62FL.csv") %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::select(c("dhsid", "dhscc", "dhsyear", "dhsclust", "rainfall_2015", "mean_temperature_2015")) %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


#---------------------------------------------------------------------------------
# write out small objects 
#---------------------------------------------------------------------------------

# write out joined HIV recode to PR, can use this for panplasmodium results
if(!dir.exists(paths = "data")){
  dir.create("data")
}

# set up figure dir for later storage
if(!dir.exists("figures")){
  dir.create("figures")
}


save(gc, file = "data/vividspace_raw.rda")
save(DRCprov, drc_stamen_back_terrain, file = "data/vividmaps_small.rda")
save(trunkroadsosm, primaryroadsosm, file = "data/osm_roads.rda")
save(riverosm, majriver, file = "data/osm_rivers.rda")
#---------------------------------------------------------------------------------
# write out large objects
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
#  geom_sf(data = majriver, color = "#9ecae1", size = 2, alpha = 0.9),
  geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true"),
  vivid_theme,
  theme(plot.background = element_blank(),
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
 #  geom_sf(data = majriver, color = "#9ecae1", size = 2, alpha = 0.9), 
 # geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true"),
  vivid_theme,
  theme(plot.background = element_blank(),
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
  #  geom_sf(data = majriver, color = "#9ecae1", size = 2, alpha = 0.9), 
  geom_sf(data = DRCprov, fill = "NA"),
  coord_sf(xlim = c(st_bbox(DRCprov)['xmin'], st_bbox(DRCprov)['xmax']), 
           ylim = c(st_bbox(DRCprov)['ymin'], st_bbox(DRCprov)['ymax']), 
           datum = NA),
  ggspatial::annotation_north_arrow(location = "bl", which_north = "true"),
  vivid_theme,
  theme(plot.background = element_blank(),
        axis.title = element_blank()) # overwrite vivid theme
)

save(prettybasemap_terraincolors, prettybasemap_hillgrey, prettybasemap_nodrc, 
     file = "data/vividmaps_large.rda")




