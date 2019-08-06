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
moranI.car.mod.framework <- spdep::moran.mc(mp$data[[2]]$plsmdprev,
                                            listw = spdep::mat2listw(car.mod.framework$W[[1]]),
                                            alternative = "greater",
                                            nsim = 1e5)

# ^ note data and W are same for all models, so can use first 



#........................
# Subset to Best Model for Plotting
#........................
bestmodel.car.mod.framework <- DIC.car.mod.framework$name[DIC.car.mod.framework$DIC == min(DIC.car.mod.framework$DIC)]
bestmodel.car.mod.framework <- car.mod.framework %>% 
  dplyr::filter(name == bestmodel.car.mod.framework)

fitted_cases_count <- as.data.frame(bestmodel.car.mod.framework$fitted)
colnames(fitted_cases_count) <- dimnames(bestmodel.car.mod.framework$W[[1]])[[1]] 

# go to long format
fitted_cases_count.ret <- fitted_cases_count %>% 
  tidyr::gather(., key = "adm1name", value = "fittedval") %>% 
  dplyr::group_by(adm1name) %>% 
  dplyr::summarise(
    pv18sprevfitted_low = quantile(fittedval, 0.025),
    pv18sprevfitted_median = quantile(fittedval, 0.5),
    pv18sprevfitted_mean = mean(fittedval),
    pv18sprevfitted_upp = quantile(fittedval, 0.975)
    
  )

#........................
# Bring it home
#........................
pv18s.adm1 <- mp$data[mp$plsmdmspec == "pv18s" & mp$maplvl == "adm1name"] 
pv18s.adm1 <- pv18s.adm1[[1]]
pv18s.adm1 <- dplyr::left_join(pv18s.adm1, fitted_cases_count.ret, by = "adm1name")

#........................
# Save out for prov
#........................
save(pv18s.adm1, car.mod.framework, 
     paramest.car.mod.framework, DIC.car.mod.framework, 
     bestmodel.car.mod.framework, 
     moranI.car.mod.framework, 
     file = "results/ProvAreal_BHM_CARBayes_models_out.rda")


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
  dplyr::filter(maplvl == "hv001" & plsmdmspec == "pv18s") %>% 
  tidyr::unnest() %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distGeo) 


pvclust <- mp$data[[5]]
sf::st_geometry(pvclust) <- NULL
pvclust.dist <- pvclust %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distGeo) 

pvclust.dist.inv <- 1/pvclust.dist
diag(pvclust.dist.inv) <- 0

mp.moranI.ret.gaus <- spdep::moran.mc(pvclust$plsmdprev,
  # temp to use MP,
                                   listw = spdep::mat2listw(pvclust.dist.inv),
                                   alternative = "greater",
  nsim = 1e5)






