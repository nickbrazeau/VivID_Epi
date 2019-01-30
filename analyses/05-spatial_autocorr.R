 #----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate spatial autocorrelations in the Pv data
# this is to determine if space is actually affecting hotspots
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(geosphere)
library(shp2graph) #Shp2graph: Tools to Convert a Spatial Network into an Igraph Graph in R -- manuscript
library(tidygraph)
library(ggraph)
library(ape)
library(spdep)
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
load("~/Documents/GitHub/VivID_Epi/data/04-basic_mapping_data.rda")
load("~/Documents/GitHub/VivID_Epi/data/osm_roads.rda")
#..........................................
# setup
#..........................................
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum", "geometry")]
clstrs <- mp %>% 
  filter(maplvl == "hv001")
clstrs$data <- map(clstrs$data, function(x){
  return( dplyr::inner_join(x, clustgeom, by = "hv001") )
})

#..........................................
# Moran's I -- several distance matrices
#..........................................
# Note, clusters are in same place for pf,pv,po so don't need to iterate over list. Can do once on any data

#.......
# greater circler
#.......
gc <- geosphere::distm(x = clstrs$data[[1]][,c("longnum", "latnum")], fun = distGeo)
gc.inv <- 1/gc
diag(gc.inv) <- 0

clstrs$MIgc <- map(clstrs$data, function(x){
  ret <- ape::Moran.I(x = x$plsmdprev, gc.inv) %>% 
    dplyr::bind_cols(.)
  return(ret)
})


#.......
# basic neighbors
#.......


#.......
# roaddistnace
#.......
rds <- sf::st_union(trunkroadsosm, primaryroadsosm)
pts <- as.matrix(clstrs$data[[1]][,c("longnum", "latnum")])[,1:2] # odd behavior bc it keep geometry ... not like a dataframe

drcroadnetwork <- shp2graph::points2network(ntdata = sf::as_Spatial(rds),
                               pointsxy = pts,
                               approach = 1) #https://rdrr.io/cran/shp2graph/man/points2network.html
# check.cnt <- nt.connect(sf::as_Spatial(rds))
# plot(check.cnt)
drcroadnetwork <- shp2graph::nel2igraph(nodelist = drcroadnetwork[[1]], edgelist = drcroadnetwork[[2]])
rds_distmat <- shortest.paths(drcroadnetwork, v=V(drcroadnetwork), to=V(drcroadnetwork))

rds_distmat.inv <- 1/rds_distmat
diag(rds_distmat.inv) <- 0

clstrs$MIrd <- map(clstrs$data, function(x){
  ret <- ape::Moran.I(x = x$plsmdprev, rds_distmat.inv) %>% 
    dplyr::bind_cols(.)
  return(ret)
})


save(drcroadnetwork, clstrs, file = "~/Documents/GitHub/VivID_Epi/data/05-spatial-autocorr.rda")







  