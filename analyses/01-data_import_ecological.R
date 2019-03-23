#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
library(sf)
library(rgeos)
library(RSAGA)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")

#---------------------------------------------------------------------------------
# set up DRC borders
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata

bb <- osmdata::getbb("Democratic Republic of the Congo", 
                     featuretype = "country",
                     format_out = "sf_polygon")

#---------------------------------------------------------------------------------
# Seasonality from Imperial/Carins
#---------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_seasonality.R")

# From OJ (via slack on 1/10/2019)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056487
# https://www.nature.com/articles/ncomms1879#supplementary-information
# https://betterexplained.com/articles/an-interactive-guide-to-the-fourier-transform/
admin_units_seasonal <- readRDS("data/imperial_share/admin1_seasonality_for_nick_from_OJ.rds")
drc_admin_units_seasonal <- admin_units_seasonal %>% 
  filter(COUNTRY_ID == "COD")


#---------------------------------------------------------------------------------
# Wetness Index
#---------------------------------------------------------------------------------
env <- rsaga.env()
env$workspace <- "~/Documents/GitHub/VivID_Epi/data/raw_data/saga_analysis/"

drcdem <- elevatr::get_elev_raster(sf::as_Spatial(bb), z=6)
raster::writeRaster(drcdem, filename = "data/raw_data/saga_analysis/drcdem.sgrd", 
                    format = "SAGA", overwrite = T)

rsaga.wetness.index(in.dem = "drcdem", 
                    out.wetness.index = "drcdem_wetnessindex", 
                    env = env)
sagawetness <- raster(x = "data/raw_data/saga_analysis/drcdem.sdat") # read in wetness

# smooth and liftover raster for ggplot 
sagawetness.pts  <-  raster::rasterToPoints(sagawetness)
sagawetness.pts.df <-  data.frame(lon = sagawetness.pts[,1], 
                                  lat = sagawetness.pts[,2], 
                                  wetness = sagawetness.pts[,3])
ggplot() + 
  geom_raster(data = sagawetness.pts.df, 
              aes(lon, lat, fill = wetness), alpha = 0.5) +
  geom_sf(data=DRCprov) +
  scale_fill_gradient2("Wetness", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") 



#---------------------------------------------------------------------------------
# Precipation Data
#---------------------------------------------------------------------------------













#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
ape <- sf::read_sf("data/redlist_species_data_primate/data_0.shp")
ape <- ape[grepl("pan paniscus|pan troglodytes|gorilla", tolower(ape$BINOMIAL)), ] %>%  # pan trog, pan panisus, gorilla sp
  dplyr::rename(species = BINOMIAL)
drc_ape <- sf::st_crop(x = ape, y = bb_sf)













save(drc_ape, file = "data/raw_data/ecological_imports.rda")