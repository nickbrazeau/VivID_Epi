#----------------------------------------------------------------------------------------------------
# Purpose of this script is to explore basic maps of Plasmodium infections in the CD2013 data
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(geosphere)
library(spdep)
library(srvyr) #wrap the survey package in dplyr syntax
library(RColorBrewer)


#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dtsrvy <- makecd2013survey(survey = dt)
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
load("data/map_bases/vivid_maps_bases.rda")

#----------------------------------------------------------------------------------------------------
# Plasmodium Point Prevalence Maps (Province & Cluster)
#----------------------------------------------------------------------------------------------------

pfldhprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pfldh) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "adm1name")
pv18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pv18s)  %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "adm1name")
po18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = po18s) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "adm1name")

pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001")
pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001")
po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001")



# bind those to a tibble
mp <- dplyr::bind_rows(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust) %>% 
  dplyr::group_by(plsmdmspec, maplvl) %>% 
  tidyr::nest()


# this awful hack becuase of this issue https://github.com/tidyverse/dplyr/issues/3483
# we are going down the rabbit hole just to try and make this stupid survey and purr package work. fine for now but return
mp$data <- lapply(list(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust), function(x) return(x))


#.............................
# Plot Summary/Point Est Maps
#..............................
pntestmaps <- pmap(mp, mapplotter)
pntestmaps <- map(pntestmaps, function(x){return(x + prettybasemap_nodrc)})



#......................
# Plot Cases
#......................
caseprevmaps <- pmap(list(data = clusters$data, plsmdmspec = clusters$plsmdmspec), casemapplotter)
caseprevmaps <- map(caseprevmaps, function(x){return(x + prettybasemap_nodrc)})



 
 #----------------------------------------------------------------------------------------------------
 # Ape Map Distributions
 #----------------------------------------------------------------------------------------------------
ape <- readRDS("data/redlist_species_data_primate/drc_ape.rds")
aperange_nhapv <- ggplot() +
   geom_sf(data = DRCprov, fill = "#d9d9d9") +
   geom_point(data = left_join(mp$data[[5]], dt[,c("hv001", "latnum", "longnum")], by = "hv001"), 
              aes(x=longnum, y=latnum, colour = plsmdprev, size = n), alpha = 0.8) +
   scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
   scale_size(guide = 'none') +
   geom_sf(data = ape, aes(fill = species), alpha = 0.4) +
   scale_fill_manual("Non-Human \n Ape Range", values = c("#33a02c", "#b3de69", "#8dd3c7", "#80b1d3")) +
   prettybasemap_nodrc +
   theme(legend.position = "bottom",
         plot.background = element_blank())





#----------------------------------------------------------------------------------------------------
# SPATIAL AUTOCORRELATION
#----------------------------------------------------------------------------------------------------

#.................................................................
# Moran's I -- greater circler
#.................................................................

gc <- geosphere::distm(x = mp$data[[1]][,c("longnum", "latnum")], fun = geosphere::distGeo) # JUST USE st_distance -- same return as geosphre
gc.inv <- 1/gc
diag(gc.inv) <- 0

mp$moranI <- lapply(mp$plsmdprev, moran.test, listw = spdep::mat2listw(gc.inv), alternative = "two.sided")
# Note, Cedar ran this in geoda for Pv and got a similar answer



 #----------------------------------------------------------------------------------------------------
 # Save Objects & Write out
 #----------------------------------------------------------------------------------------------------
 
saveRDS(mp, file = "data/derived_data/basic_cluster_mapping_data.rds")

save(pntestmaps, caseprevmaps, aperange_nhapv,
      file = "results/basic_maps_results.rda")
