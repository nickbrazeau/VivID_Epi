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

#----------------------------------------------------------------------------------------------------
# Plasmodium Point Prevalence Maps (Province & Cluster)
#----------------------------------------------------------------------------------------------------

pfldhprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pfldh, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "adm1name")
pv18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pv18s, adm1shp = DRCprov)  %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "adm1name")
po18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = po18s, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "adm1name")

pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001")
pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001")
po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s, adm1shp = DRCprov) %>% 
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


#......................
# Plot Cases
#......................
clst <- mp %>% 
  dplyr::filter(maplvl == "hv001") %>% 
  dplyr::select(-c("maplvl"))

caseprevmaps <- purrr::pmap(list(data = clst$data, plsmdmspec = clst$plsmdmspec), casemap_prev_plotter)

case_n_maps <- purrr::pmap(list(data = clst$data, plsmdmspec = clst$plsmdmspec), casemap_n_plotter)




#----------------------------------------------------------------------------------------------------
# Ape Map Distribution for Pv
#----------------------------------------------------------------------------------------------------
ape <- readRDS("data/derived_data/drc_ape.rds")
aperange_nhapv <- 
  caseprevmaps[[2]] +
  geom_sf(data = ape, aes(fill = species), alpha = 0.4) +
  scale_fill_manual("Non-Human \n Ape Range", values = c("#33a02c", "#b3de69", "#8dd3c7", "#80b1d3"))





 #----------------------------------------------------------------------------------------------------
 # Save Objects & Write out
 #----------------------------------------------------------------------------------------------------
 
saveRDS(mp, file = "data/derived_data/basic_cluster_mapping_data.rds")

save(pntestmaps, caseprevmaps, case_n_maps, aperange_nhapv,
      file = "results/basic_maps_results.rda")
