#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make ICAR models
# Will aggregate up to the province level
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/R/00-functions_epi.R") 
library(tidyverse)
library(srvyr) #wrap the survey package in dplyr syntax
library(CARBayes)

#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dtsrvy <- makecd2013survey(survey = dt)
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")


#......................
# Subset to Pv
#......................
pvprov <- mp %>% 
  dplyr::filter(plsmdmspec == "pv18s" & maplvl == "adm1name") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pvprov <- sf::st_as_sf(pvprov)

#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------

prov.covar <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(adm1name) %>% 
  dplyr::summarise(n = srvyr::survey_total(count), 
                   meanprov_precip_lag_cont_scale_clst = srvyr::survey_mean(precip_lag_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95),
                   meanprov_alt_dem_cont_scale_clst = srvyr::survey_mean(alt_dem_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95)
  )

pvprov <- dplyr::left_join(pvprov, prov.covar, by = "adm1name")

#-------------------------------------------------------------------------
# Basic Multilevel Model
#-------------------------------------------------------------------------
# https://github.com/tmalsburg/MCMCglmm-intro
library(MCMCglmm)



prior <- list(
  R=list(V=1, n=1, fix=1),
  G=list(G1=list(V        = diag(8),
                 n        = 8,
                 alpha.mu = rep(0, 8),
                 alpha.V  = diag(8)*25^2),
         G2=list(V        = diag(4),
                 n        = 4,
                 alpha.mu = rep(0, 4),
                 alpha.V  = diag(4)*25^2)))

m3 <- MCMCglmm(pronoun ~ (a + b + c)^3,
               ~ us(1 + (a + b + c)^3):subject +
                 us(1 + (a + b    )^2):item,
               data   = d,
               family = "categorical",
               prior  = prior.m3,
               thin   = 1,
               burnin = 3000,
               nitt   = 4000)


#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model 
#-------------------------------------------------------------------------
#......................
# Make Adjacency Matrix for Pv 
#......................
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
W.nb <- spdep::poly2nb(sf::as_Spatial(pvprov), row.names = rownames(sf::as_Spatial(pvprov)@data), style = "B") # binary weights taking values zero or one (only one is recorded)












