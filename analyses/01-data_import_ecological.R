#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
devtools::install_github("malaria-atlas-project/malariaAtlas")
library(malariaAtlas)
devtools::install_github("OJWatson/magenta")
library(magenta)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")

#---------------------------------------------------------------------------------
# set up DRC borders
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata

bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')

#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
primate <- sf::read_sf("data/redlist_species_data_primate/data_0.shp")
ape <- primate[grepl("pan paniscus|pan troglodytes|gorilla", tolower(primate$BINOMIAL)), ] # pan trog, pan panisus, gorilla sp


#---------------------------------------------------------------------------------
# Seasonality from Imperial/Carins
#---------------------------------------------------------------------------------

# From OJ (via slack on 1/10/2019)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056487
# https://www.nature.com/articles/ncomms1879#supplementary-information
# https://betterexplained.com/articles/an-interactive-guide-to-the-fourier-transform/
admin_units_seasonal <- readRDS("data/imperial_share/admin1_seasonality_for_nick_from_OJ.rds")
drc_admin_units_seasonal <- admin_units_seasonal %>% 
  filter(COUNTRY_ID == "COD")

### NEED TO FINISH -- ERROR WITH THETA_C still
source("~/Documents/GitHub/VivID_Epi/R/00-functions_seasonality.R")




#---------------------------------------------------------------------------------
# Rainfall Data from OJ for admin1 in DHS
#---------------------------------------------------------------------------------
# Daily DRC's real rainfall data from CHIRPs for 2013-2017

rain <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/imperial_share/drc_rainfall_2013-2017.rds")

 
#---------------------------------------------------------------------------------
# MAP project for R 
#---------------------------------------------------------------------------------
# from their example code chunk
# download raster of travel time to cities (Weiss et al 2018) for study area & visualise this
TravelToCities <- malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                                  extent = bb)
TravelToCities <- log(TravelToCities+0.1) # MAP project decided 0.1 for offest. going to keep it



#---------------------------------------------------------------------------------
# Link to GE
#---------------------------------------------------------------------------------








