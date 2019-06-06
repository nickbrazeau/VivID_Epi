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
# Precipation Data
#---------------------------------------------------------------------------------
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
ape <- sf::read_sf("data/redlist_species_data_primate/data_0.shp")
ape <- ape[grepl("pan paniscus|pan troglodytes|gorilla", tolower(ape$BINOMIAL)), ] %>%  # pan trog, pan panisus, gorilla sp
  dplyr::rename(species = BINOMIAL)
drc_ape <- sf::st_crop(x = ape, y = bb)


saveRDS(drc_ape, file = "data/redlist_species_data_primate/drc_ape.rds")
