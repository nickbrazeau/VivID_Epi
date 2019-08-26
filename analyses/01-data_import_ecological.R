#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(sf)



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


saveRDS(drc_ape, file = "data/derived_data/drc_ape.rds")




#---------------------------------------------------------------------------------
# Precipation Data
#---------------------------------------------------------------------------------
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

# note, no hiv testing in 70; 271/318 lost to contamination


#----------------------------------------------------------------------------------------------------
# Waterways
#----------------------------------------------------------------------------------------------------
wtrlns <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_lines_shp/hotosm_cod_waterways_lines.shp") %>% 
  dplyr::filter(waterway %in% c("stream", "river", "riverbank")) %>% 
  dplyr::rename(watertype = waterway) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))

wtrply <- sf::read_sf("data/raw_data/hotosm_data/hotosm_cod_waterways_polygons_shp/hotosm_cod_waterways_polygons.shp") %>% 
  dplyr::filter(water == "lake") %>% 
  dplyr::rename(watertype = water) %>% 
  dplyr::select(c("osm_id", "watertype", "geometry"))

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
hlthsites.harvard.drc <- readxl::read_excel("data/raw_data/harvard_dataverse/Ouma_Okiro_Snow_Africa_Hospitals_Data.xlsx") %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  magrittr::set_colnames(gsub(pattern = " ", "_", colnames(.))) %>% 
  dplyr::filter(country == "Democratic Republic of Congo") %>%
  sf::st_as_sf(coords = c("long", "lat"), 
               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))


#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

ge.osrm <- ge %>% 
  dplyr::select("dhsclust") 
rownames(ge.osrm) <- ge.osrm$dhsclust


hlthsites.harvard.drc.osrm <- hlthsites.harvard.drc %>% 
  dplyr::mutate(id = paste0("hlth", seq(1:nrow(.))),
                id = factor(id)) %>%
  dplyr::select(id)

rownames(hlthsites.harvard.drc.osrm) <- hlthsites.harvard.drc.osrm$id


# now access API
hosptrvltmes <- osrm::osrmTable(src = ge.osrm,
                                dst = hlthsites.harvard.drc.osrm,
                                measure = "duration") # in meters


hlthdist_out <- 
  cbind.data.frame(hv001 = as.numeric( rownames(hosptrvltmes$duration) ), hosptrvltmes$duration ) %>% 
  tidyr::gather(., key = "hsptl", value = "hlthst_duration", 2:ncol(.)) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    hlthst_nrst_duration = min(hlthst_duration)
  ) 



#..........
# Note, cluter 469 cannot be resolved with osrm
# however cluster 313 is nearby, approximately 19763.25 
# just going to add this greater circle distance
#..........
clst469 <- ge %>% 
  dplyr::filter(dhsclust == 469)

dist <- raster::pointDistance(p1 = sf::as_Spatial(clst469),
                              p2 = sf::as_Spatial(ge),
                              lonlat = T)
# find nearest cluster
nghbrfor469 <- ge$dhsclust[which(dist %in%  sort(dist)[2])]
hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469] 
sort(dist)[2] # need to travel about 20 km, which is about 1hr by car in rural setting

hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == 469] <- hlthdist_out$hlthst_nrst_duration[hlthdist_out$hv001 == nghbrfor469] + 60


#----------------------------------------------------------------------------------------------------
# write out
#----------------------------------------------------------------------------------------------------
saveRDS(object = wtrdist_out, file = "data/derived_data/hotosm_waterways_dist.rds")
saveRDS(object = hlthdist_out, file = "data/derived_data/hlthdist_out_minduration.rds")
save(wtrply, wtrlns, file = "data/derived_data/hotosm_waterways.RDA")
