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





##################################################################################
##########                    Impute Precipitation                ################
##################################################################################
# find missing clusters
clst.all <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds") %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::filter(hv103 == 1) %>% 
  dplyr::select(c("hv001", "annual_precipitation_2015", "latnum", "longnum", "geometry")) %>% 
  dplyr::filter(!duplicated(.))

missclust <- clst.all[clst.all$annual_precipitation_2015 == -9999, ]
notmissclust <- clst.all[clst.all$annual_precipitation_2015 != -9999, ]
# can we do k-nearest neighbors average again
ggplot() + 
  geom_sf(data=notmissclust, color= "black") + 
  geom_sf(data=missclust, color= "red", size = 2, alpha  = 0.5) 

# first missing cluster
dist.clust1 <- raster::pointDistance(p1 = sf::as_Spatial(missclust[1,]),
                              p2 = sf::as_Spatial(notmissclust),
                              lonlat = T)
# find 5 nearby clusters
dist.clust1.sorted.5 <- sort(dist.clust1)[1:5]
nrbyclstrs.clust1 <- which(dist.clust1 %in% dist.clust1.sorted.5)


# SECOND missing cluster
dist.clust2 <- raster::pointDistance(p1 = sf::as_Spatial(missclust[2,]),
                                     p2 = sf::as_Spatial(notmissclust),
                                     lonlat = T)
# find 5 nearby clusters
dist.clust2.sorted.5 <- sort(dist.clust2)[1:5]
nrbyclstrs.clust2 <- which(dist.clust2 %in% dist.clust2.sorted.5)

# mean precip
precip1 <- mean( notmissclust$annual_precipitation_2015[notmissclust$hv001 %in% nrbyclstrs.clust1] )
precip2 <- mean( notmissclust$annual_precipitation_2015[notmissclust$hv001 %in% nrbyclstrs.clust2] )

# overwrite
clst.all$annual_precipitation_2015[clst.all$annual_precipitation_2015 == -9999] <- c(precip1, precip2)

# rename
clst.all <- clst.all %>% 
  dplyr::rename(annual_precipitation_2015imp = annual_precipitation_2015) %>% 
  dplyr::select(c("hv001", "annual_precipitation_2015imp"))
sf::st_geometry(clst.all) <- NULL

saveRDS(object = clst.all, file = "data/derived_data/annual_precipitation_2015_imputed.RDS")


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
