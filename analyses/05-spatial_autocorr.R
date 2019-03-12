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
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R")
 
 
 #......................
 # Import Data
 #......................
 dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
 options(survey.lonely.psu="certainty")
 dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, strata = hv023, weights = hv005_wi)
 
 #spatial from the DHS -- these are cluster level vars
 ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
 colnames(ge) <- tolower(colnames(ge))
 ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
 ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
 # remove clusters that were missing from the DHS, see readme
 ge <- ge %>% 
   dplyr::rename(hv001 = dhsclust) # for easier merge with PR
 
 #......................
 # Summarize by Cluster
 #......................
 pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh) %>% 
   dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001") %>% 
   dplyr::left_join(x=., y = ge)
 pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s) %>% 
   dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001") %>% 
   dplyr::left_join(x=., y = ge)
 po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s) %>% 
   dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001") %>% 
   dplyr::left_join(x=., y = ge)
 
#..........................................
# Moran's I -- several distance matrices
#..........................................
# Note, clusters are in same place for pf,pv,po so don't need to iterate over list. Can do once on any data
 clstrs <- dplyr::bind_rows(pfldhclust, pv18sclust, po18sclust) %>% 
   dplyr::group_by(plsmdmspec) %>% 
   tidyr::nest()
 
 # this awful hack becuase of this issue https://github.com/tidyverse/dplyr/issues/3483
 # we are going down the rabbit hole just to try and make this stupid survey and purr package work. fine for now but return
 clstrs$data <- lapply(list(pfldhclust, pv18sclust, po18sclust), function(x) return(x))
 
 
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

moran.test(clstrs$data[[2]]$plsmdprev, mat2listw(gc.inv), 
           alternative = "two.sided")


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







  