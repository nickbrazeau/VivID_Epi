#..............................................................
# Purpose of this script is to have a backend for the CarBayes models
#..............................................................
source("R/00-functions_basic.R") 
source("R/00-functions_epi.R") 
source("R/00-MCMC_diagnostics.R")
library(tidyverse)
library(furrr)
library(srvyr) #wrap the survey package in dplyr syntax
library(CARBayes)
library(HDInterval)
library(raster)
set.seed(48)

#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")

#------------------------------------------------------------------------
# Subset to Pv
#------------------------------------------------------------------------
pvprov.weighted.nosf <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(adm1name) %>% 
  dplyr::summarise(n = srvyr::survey_total(count), 
                   plsmdn = srvyr::survey_total(pv18s, na.rm = T), 
                   plsmdprev = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95))

# need to keep integers, so will round
pvprov.weighted.nosf <- pvprov.weighted.nosf %>% 
  dplyr::mutate(plsmdn = round(plsmdn, 0),
                n = round(n, 0))



#..............................................................
# Import the Covariates
#..............................................................
pvcovar <- readRDS("data/derived_data/vividepi_prov_covars_bayesian_fit.RDS")
# combine
pvprov.weighted.nosf <- dplyr::left_join(pvprov.weighted.nosf, pvcovar, by = "adm1name")



############################################################################################################
#######################                      Diagnostic Chains                       #######################             
#############################################################################################################

#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model 
#-------------------------------------------------------------------------
#......................
# Make Adjacency Matrix for Pv 
#......................
DRCprovsp <- dplyr::left_join(pvprov.weighted.nosf, DRCprov)
DRCprovsp <- sf::st_as_sf(DRCprovsp)
# sanity
DRCprovsp <- sf::st_transform(DRCprovsp, crs = "+init=epsg:4326")
sp::identicalCRS(sf::as_Spatial(DRCprovsp), sf::as_Spatial(DRCprov))

# make adj
W.nb <- spdep::poly2nb(sf::as_Spatial(DRCprovsp), row.names = DRCprovsp$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

#......................
# Make Model Framework
#......................
prov.covar.names <- c("precip_scale", "crop_scale", "hlthdist_scale")
mod.framework <- tibble(name = c("CAR_intercept", "CAR_covar"),
                        formula = c("plsmdn ~ 1", 
                                    paste0("plsmdn ~ ", paste(prov.covar.names, collapse = " + "))),
                        burnin = 1e4,
                        n.sample = 1e4+1e4,
                        family = "binomial"
)


mod.framework$trials <- lapply(1:nrow(mod.framework), function(x) return(pvprov.weighted.nosf$n))
mod.framework$data <- lapply(1:nrow(mod.framework), function(x) return(pvprov.weighted.nosf))
mod.framework$W <- lapply(1:nrow(mod.framework), function(x) return(W))
mod.framework$thin <- 1 # no thinning for diagnostics



# replicate this four times for our four chains
mod.framework <- lapply(1:4, function(x) return(mod.framework)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(name)

#......................
# Make a wrapper for CARBAYES
#......................
wrap_S.CARleroux <- function(name, formula, family, trials, W, data, burnin, n.sample, thin){
  
  formvec <- paste(formula, collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 2 # need intercept and betas
  betacount <- ifelse(grepl("1", formvec), 1, betacount) # corner case of just intercept
  prior.var.betavec <- rep(100, betacount) # note prior default of 1e5 setting here
  
  # rho here is NULL by default and esimtated in the model
  ret <- CARBayes::S.CARleroux(formula = as.formula(formula), 
                               family = family, 
                               trials = trials, 
                               W = W,
                               data = data,
                               burnin = burnin, 
                               prior.var.beta = prior.var.betavec,
                               n.sample = n.sample,
                               thin = thin)
  
  
  return(ret)
}


mod.framework$MCMC <- purrr::pmap(mod.framework, wrap_S.CARleroux)

#..............................................................
# save out diagnostic chains
#..............................................................
dir.create("analyses/07-spatial_prediction/ProvModels/", recursive = T)
saveRDS(mod.framework, "analyses/07-spatial_prediction/ProvModels/ProvModel_diag_chains.RDS")


############################################################################################################
#######################                         Long Chains                          #######################             
#############################################################################################################
mod.framework.long <- mod.framework[,c("name", "formula", "family", "trials", "W", "data", "burnin", "n.sample")] %>% 
  dplyr::filter(!duplicated(.))

mod.framework.long$burnin <- 1e5
mod.framework.long$n.sample <- 1e5 + 1e5
mod.framework.long$thin <- 10 # some thing for longer run 

mod.framework.long$MCMC <- purrr::pmap(mod.framework.long, wrap_S.CARleroux)

#..............................................................
# save out long chains
#..............................................................
saveRDS(mod.framework.long, "analyses/07-spatial_prediction/ProvModels/ProvModel_long_chains.RDS")




