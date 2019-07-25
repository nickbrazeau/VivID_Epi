#----------------------------------------------------------------------------------------------------
# Purpose of this script is to make ICAR models aggregated at the Province Level
# We will run four chains for each to see if the are all converging in the same place
# and appear to be mixing well
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
W.nb <- spdep::poly2nb(sf::as_Spatial(pvprov.weighted), row.names = pvprov.weighted$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

#......................
# Make Model Framework
#......................
prov.covar.names <- c("meanprov_precip_lag_cont_scale_clst", "meanprov_alt_dem_cont_scale_clst")
mod.framework <- tibble(name = c("riid", "ICAR", "CAR", "riid_covar", "ICAR_covar", "CAR_covar"),
                        formula = c("plsmdn ~ 1", "plsmdn ~ 1", "plsmdn ~ 1",
                                    rep(paste0("plsmdn ~ ", paste(prov.covar.names, collapse = " + ")), 3)),
                        rho = c(0, 1, NA, 0, 1, NA),
                        burnin = 1e3,
                        n.sample = 1e5,
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
wrap_S.CARleroux <- function(formula, family, trials, W, rho, data, burnin, n.sample){
  if(!is.na(rho)){
    ret <- CARBayes::S.CARleroux(formula = as.formula(formula), 
                                 family = family, 
                                 trials = trials, 
                                 W = W,
                                 rho = rho,
                                 data = data,
                                 burnin = burnin, 
                                 n.sample = n.sample)
  } else if(is.na(rho)){ # rho needs to be NULL which has a hard time in a vector
    ret <- CARBayes::S.CARleroux(formula = as.formula(formula), 
                                 family = family, 
                                 trials = trials, 
                                 W = W,
                                 data = data,
                                 burnin = burnin, 
                                 n.sample = n.sample)
  }
  
  return(ret)
}


mod.framework$MCMC <- purrr::pmap(mod.framework[, -1], wrap_S.CARleroux)



#-------------------------------------------------------------------------
# MCMC Diagnostics
#-------------------------------------------------------------------------
mod.framework$mcmc.modsum <- purrr::map(mod.framework$MCMC, print) # note, print is overloaded here
mod.framework$summresults <- purrr::map(mod.framework$mcmc.modsum, "summary.results")
mod.framework$summresults[1:4]
mod.framework$summresults[5:8]
mod.framework$summresults[9:12]
mod.framework$summresults[13:16]
mod.framework$summresults[17:20]
mod.framework$summresults[21:24]


#............................
# Let's Look at the Chains
#............................
mod.framework$samples <- purrr::map(mod.framework$MCMC, "samples")
mod.framework$beta <- purrr::map(mod.framework$samples, "beta")
mod.framework$phi <- purrr::map(mod.framework$samples, "phi")
mod.framework$tau2 <- purrr::map(mod.framework$samples, "tau2")
mod.framework$rho <- purrr::map(mod.framework$samples, "rho")
mod.framework$fitted <- purrr::map(mod.framework$samples, "fitted")


chains <- mod.framework %>% 
  dplyr::select(c("name", "beta", "phi", "tau2", "rho", "fitted")) %>% 
  dplyr::group_by(name) %>% 
  tidyr::nest()


make_mcmc_chain_plots <- function(chaindat, tempdir){
  jpeg(paste0(tempdir, "_", names(chaindat), ".jpeg"), height = 8, width = 11, units = "in", res=200)
  par(mfrow=c(2,2))
  plot(chaindat[[1]])
  plot(chaindat[[2]])
  plot(chaindat[[3]])
  plot(chaindat[[4]])
  graphics.off()

}
# mytempdir <- "~/Desktop/nfbtemp/"
# purrr::map(chains$data, function(x){
#  apply(x, 2, make_mcmc_chain_plots, tempdir = mytempdir)
# })


#-------------------------------------------------------------------------
# Take Very Long Chains to LL 
#-------------------------------------------------------------------------
setwd("analyses/08-spatial_prediction/areal_results/")
mod.framework.slurm <- mod.framework[,c("formula", "family", "trials", "W", "rho", "data", "burnin", "n.sample")] %>% 
  dplyr::filter(!duplicated(.))

mod.framework.slurm$burnin <- 1e5
mod.framework.slurm$n.sample <- 1e8

ntry <- nrow(mod.framework.slurm)
sjob <- rslurm::slurm_apply(f = wrap_S.CARleroux, 
                            params = mod.framework.slurm, 
                            jobname = 'MCMCCAR_models',
                            nodes = 1, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 128000,
                                                 array = sprintf("0-", 
                                                                 ntry),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5-00:00:00"))

cat("*************************** \n Submitted CARBayes Models \n *************************** ")







