 #----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate spatial autocorrelations in the Pv data
# this is to determine if space is actually affecting hotspots
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions_basic.R") 
library(tidyverse)
library(sf)
library(geosphere)
library(ape)
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
load("~/Documents/GitHub/VivID_Epi/data/04-basic_mapping_data.rda")
#..........................................
# setup
#..........................................
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum")]
clusters <- mp %>% 
  filter(maplvl == "hv001")
clusters$data <- map(clusters$data, function(x){
  return( dplyr::inner_join(x, clustgeom, by = "hv001") )
})

#..........................................
# Moran's I -- several distance matrices
#..........................................
# Note, clusters are in same place for pf,pv,po so don't need to iterate over list. Can do once on any data
clusters$gc <- geosphere::distHaversine(clusters$data[[1]][,c("longnum", "latnum")])





  