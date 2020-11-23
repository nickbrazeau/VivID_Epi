#..............................................................
# Purpose of this script is to have a backend for the CarBayes models
#..............................................................
source("R/00-functions_basic.R") 
source("R/00-functions_epi.R") 
source("R/00-MCMC_diagnostics.R")
library(tidyverse)
library(srvyr) #wrap the survey package in dplyr syntax
library(CARBayes)
library(HDInterval)
library(raster)
set.seed(48)
tol <- 1e-3

#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")


#------------------------------------------------------------------------
# Subset to Pv
#------------------------------------------------------------------------
pvprov.weighted <- mp$data[mp$plsmdmspec == "pv18s" & mp$maplvl == "adm1name"][[1]]
# vectors have destroyed spatial class, need to remake
pvprov.weighted <- sf::st_as_sf(pvprov.weighted)
sf::st_crs(pvprov.weighted) <-  sf::st_crs(ge)
# need ints (binomail prob), so will round
pvprov.weighted <- pvprov.weighted %>% 
  dplyr::mutate(plsmdn = round(plsmdn, 0),
                n = round(n, 0))

pvprov.weighted.nosf <- pvprov.weighted
sf::st_geometry(pvprov.weighted.nosf) <- NULL

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
W.nb <- spdep::poly2nb(sf::as_Spatial(pvprov.weighted), row.names = pvprov.weighted$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

#......................
# Make Model Framework
#......................
prov.covar.names <- c("precip_scale", "crop_scale", "nightlight_scale")
mod.framework <- tibble(name = c("CAR_intercept", "CAR_covar"),
                        formula = c("plsmdn ~ 1", 
                                    paste0("plsmdn ~ ", paste(prov.covar.names, collapse = " + "))),
                        burnin = 1e3,
                        n.sample = 1e4 + 1e3,
                        family = "binomial"
)


mod.framework$trials <- lapply(1:nrow(mod.framework), function(x) return(pvprov.weighted.nosf$n))
mod.framework$data <- lapply(1:nrow(mod.framework), function(x) return(pvprov.weighted.nosf))
mod.framework$W <- lapply(1:nrow(mod.framework), function(x) return(W))



# replicate this four times for our four chains
mod.framework <- lapply(1:4, function(x) return(mod.framework)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(name)

#......................
# Make a wrapper for CARBAYES
#......................
wrap_S.CARleroux <- function(name, formula, family, trials, W, data, burnin, n.sample){
  
  formvec <- paste(formula, collapse = "")
  betacount <- stringr::str_count(formvec, "\\+") + 2 # need intercept and betas
  betacount <- ifelse(grepl("1", formvec), 1, betacount) # corner case of just intercept
  prior.var.betavec <- rep(5e4, betacount) # note prior setting here
  
  # rho here is NULL by default and esimtated in the model
  ret <- CARBayes::S.CARleroux(formula = as.formula(formula), 
                               family = family, 
                               trials = trials, 
                               W = W,
                               data = data,
                               burnin = burnin, 
                               prior.var.beta = prior.var.betavec,
                               n.sample = n.sample)
  
  
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

mod.framework.long$burnin <- 1e3
mod.framework.long$n.sample <- 1e5 + 1e3

mod.framework.long$MCMC <- purrr::pmap(mod.framework.long, wrap_S.CARleroux)

#..............................................................
# save out long chains
#..............................................................
saveRDS(mod.framework.long, "analyses/07-spatial_prediction/ProvModels/ProvModel_long_chains.RDS")




