#----------------------------------------------------------------------------------------------------
# Purpose of this script is to collect and plot the fitted/predicted values from the CAR model
# and gaussian process model 
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R") 
source("R/00-functions_epi.R") 
source("R/00-MCMC_diagnostics.R")
library(tidyverse)
library(CARBayes)
library(PrevMap)
library(sf)
load("data/map_bases/vivid_maps_bases.rda")
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")

##########################################################################################################################################
##########################################################################################################################################
##############                     Province Model                      ###################################################################                     
##########################################################################################################################################
##########################################################################################################################################
# readRDS
car.mod.framework$mcmc.modsum <- purrr::map(car.mod.framework$MCMC, print) 
car.mod.framework$summresults <- purrr::map(car.mod.framework$mcmc.modsum, "summary.results")
car.mod.framework$modfit <- purrr::map(car.mod.framework$mcmc.modsum, "modelfit")
car.mod.framework$DIC <- purrr::map(car.mod.framework$modfit, "DIC")

#........................
# Check CARBayes Chains
#........................
car.mod.framework$samples <- purrr::map(car.mod.framework$MCMC, "samples")
car.mod.framework$beta <- purrr::map(car.mod.framework$samples, "beta")
car.mod.framework$phi <- purrr::map(car.mod.framework$samples, "phi")
car.mod.framework$tau2 <- purrr::map(car.mod.framework$samples, "tau2")
car.mod.framework$rho <- purrr::map(car.mod.framework$samples, "rho")
car.mod.framework$fitted <- purrr::map(car.mod.framework$samples, "fitted")


chains <- car.mod.framework %>% 
  dplyr::select(c("name", "beta", "phi", "tau2", "rho", "fitted")) %>% 
  dplyr::group_by(name) %>% 
  tidyr::nest()


mytempdir <- "analyses/08-spatial_prediction/prov_map_long_chains/"
dir.create(mytempdir, recursive = T)
wrap_chain_plotter(tempdir = mytempdir, chains = chains)

#........................
# Parameter Est Summaries
#........................
names(car.mod.framework$summresults) <- car.mod.framework$name
car.mod.framework$summresults <- purrr::map(car.mod.framework$summresults, tibble::as_tibble)
paramest.car.mod.framework <- car.mod.framework$summresults %>% 
  dplyr::bind_rows(., .id = "model_formulation")


DIC.car.mod.framework <- car.mod.framework %>% 
  dplyr::select(c("name", "DIC")) %>% 
  dplyr::filter(!duplicated(.)) %>% 
  tidyr::unnest(DIC) %>% 
  dplyr::mutate(DIC = round(DIC, 3))




#........................
# Moran's I
#........................
mp.moranI.ret <- spdep::moran.mc(mp$data[[2]]$plsmdprev,
                                 listw = spdep::mat2listw(mod.framework$W[[1]]),
                                 alternative = "greater",
                                 nsim = 1e5)

# ^ note data and W are same for all models, so can use first 



#........................
# Subset to Best Model for Plotting
#........................

##########################################################################################################################################
##########################################################################################################################################
##############                     Cluster Model                      ###################################################################                     
##########################################################################################################################################
##########################################################################################################################################


# readRDS
# 

#........................
# Moran's I
#........................
gc <- mp %>% 
  dplyr::filter(maplvl == "hv001" & plsmdmspec == "pfldh") %>% # pfldh just so we can have one data obj at correct map level
  tidyr::unnest() %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distGeo) 





pvclust <- mp$data[[5]]
sf::st_geometry(pvclust) <- NULL
pvclust.dist <- pvclust %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distGeo) 

pvclust.dist <- log(pvclust.dist)
pvclust.dist.inv <- 1/pvclust.dist
diag(pvclust.dist.inv) <- 0

mp.moranI.ret.gaus <- spdep::moran.mc(pvclust$plsmdprev,
  # temp to use MP,
                                   listw = spdep::mat2listw(pvclust.dist.inv),
                                   alternative = "greater",
  nsim = 1e5)






