#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make ICAR models
# Will aggregate up to the province level
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R") 
source("R/00-functions_epi.R") 
source("R/00-MCMC_diagnostics.R")
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
pvprov.weighted <- mp %>% 
  dplyr::filter(plsmdmspec == "pv18s" & maplvl == "adm1name") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pvprov.weighted <- sf::st_as_sf(pvprov.weighted)
# need to keep integers, so will round
pvprov.weighted <- pvprov.weighted %>% 
  dplyr::mutate_if(is.numeric, round, 0)




#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------

prov.covar <- dtsrvy %>% 
  dplyr::group_by(adm1name) %>% 
  dplyr::summarise(meanprov_precip_lag_cont_scale_clst = srvyr::survey_mean(precip_lag_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95),
                   meanprov_alt_dem_cont_scale_clst = srvyr::survey_mean(alt_dem_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95)
  )

pvprov.weighted <- dplyr::left_join(pvprov.weighted, prov.covar, by = "adm1name")
pvprov.weighted.nosf <- pvprov.weighted
sf::st_geometry(pvprov.weighted.nosf) <- NULL


#-------------------------------------------------------------------------
# Basic Multilevel Model
#-------------------------------------------------------------------------
# # https://github.com/tmalsburg/MCMCglmm-intro
# # https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
# # https://ms.mcmaster.ca/~bolker/R/misc/foxchapter/bolker_chap.html
# library(MCMCglmm)
# 
# # Random Intercept model, no covariates
# prior.randint <- list(R = list(V = diag(1), nu = 1e-3),
#                       G = list(G1 = list(V = diag(1), nu = 1, 
#                                          alpha.mu = 0, alpha.V = diag(1)*1e2)))
# post.randint <- MCMCglmm(cbind(plsmdn, n) ~ 1, 
#                          random = ~adm1name,
#                          prior = prior.randint, 
#                          data = pvprov.weighted.nosf, 
#                          family = "multinomial2",
#                          nitt = 1e5,
#                          burnin = 1e3)
# 
# # check model fit
# coda::autocorr(post.randint$Sol)
# coda::effectiveSize(post.randint$Sol)
# plot(post.randint$Sol)
# 
# # pretty miserable looking, likely because too little data?
# # Within an MCMCglmm fit the chains are stored in two separate matrices (mcmc objects, actually, but these can be treated a lot like matrices) called  Sol (fixed effects) and VCV (variances and covariances). (Random effects are not saved unless you set pr=TRUE.)
# allChains <- as.mcmc(cbind(post.randint$Sol,post.randint$VCV))
# plot(allChains)

#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model 
#-------------------------------------------------------------------------
#......................
# Make Adjacency Matrix for Pv 
#......................
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
W.nb <- spdep::poly2nb(sf::as_Spatial(pvprov.weighted), row.names = pvprov$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)


#......................
# Fit Random Intercept (set rho to 0) 
#......................
rand.int <- CARBayes::S.CARleroux(formula = plsmdn~1 , 
                                  family = "binomial", 
                                  trials = pvprov.weighted.nosf$n, 
                                  W = W,
                                  rho = 0,
                                  data = pvprov.weighted.nosf,
                                  burnin = 10, 
                                  n.sample = 50)


#......................
# Fit Complete Spatial Dependence (set rho to 1) 
#......................
spat.int <- CARBayes::S.CARleroux(formula = plsmdn~1 , 
                                  family = "binomial", 
                                  trials = pvprov.weighted.nosf$n, 
                                  W = W,
                                  rho = 1,
                                  data = pvprov.weighted.nosf,
                                  burnin = 10, 
                                  n.sample = 50)



#......................
# Fit CARleroux with rho being estimated
#......................
CARleroux <- CARBayes::S.CARleroux(formula = plsmdn~1 , 
                                  family = "binomial", 
                                  trials = pvprov.weighted.nosf$n, 
                                  W = W,
                                  rho = NULL,
                                  data = pvprov.weighted.nosf,
                                  burnin = 10, 
                                  n.sample = 50)

model.spatial$fitted.values
length(model.spatial$fitted.values)
length(Y)

plot(Y ~ model.spatial$fitted.values)











