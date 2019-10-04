
source("R/00-functions_basic.R")
source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) 
library(rgeos)
library(rgdal)
library(geoR)
library(raster)
library(PrevMap)
set.seed(48)
tol <- 1e-3
#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
ge <- sf::st_as_sf( readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds") )
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")

riskvars = c("precip_mean_cont_scale_clst")

#......................
# Subset to Pv
#......................
pvclust.weighted <- mp$data[mp$plsmdmspec == "pv18s" & mp$maplvl == "hv001"][[1]]
# vectors have destroyed spatial class, need to remake
pvclust.weighted <- sf::st_as_sf(pvclust.weighted)
sf::st_crs(pvclust.weighted) <-  sf::st_crs(ge)

# need to keep integers, so will round
pvclust.weighted <- pvclust.weighted %>% 
  dplyr::mutate(plsmdn = round(plsmdn, 0),
                n = round(n, 0))



#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------

pvclst.covar <- dtsrvy %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(precip_mean_cont_scale_clst = srvyr::survey_mean(precip_mean_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95) # identical by clst but ok
  )

pvclust.weighted <- dplyr::left_join(pvclust.weighted, pvclst.covar, by = "hv001")
pvclust.weighted.nosf <- pvclust.weighted
sf::st_geometry(pvclust.weighted.nosf) <- NULL

# transform count of "successes" to logit space
pvclust.weighted.nosf$plsmdlogit <- log( (pvclust.weighted.nosf$plsmdn + 0.5)/(pvclust.weighted.nosf$n - pvclust.weighted.nosf$plsmdn + 0.5) ) # 0.5 as tolerance for 0s

#-------------------------------------------------------------------------
# Make Model Framework
#-------------------------------------------------------------------------
mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                formula = c("plsmdlogit ~ 1", 
                                            paste0("plsmdlogit ~ ", paste(riskvars, collapse = "+")))
)

coords <- as.formula(paste0("~ longnum + latnum"))

mod.framework$data <- lapply(1:nrow(mod.framework), function(x) return(pvclust.weighted.nosf))

mod.framework$coords <- lapply(1:nrow(mod.framework), function(x) return(coords))


####################################################################################
###########           PRIORS  & Direction Diagnostic                      ##########
#####################################################################################
fit.glm <- glm(cbind(plsmdn, n - plsmdn) ~ 1, 
               data = pvclust.weighted.nosf,
               family = binomial)


mypriors.intercept <- PrevMap::control.prior(beta.mean = 0,
                                             beta.covar = 1e4,
                                             log.normal.nugget = c(-5,5), # this is tau2
                                             uniform.phi = c(0,10),
                                             log.normal.sigma = c(0,25)
)

# NB covar matrix
covarsmat <- matrix(0, ncol = 2, nrow=2)
diag(covarsmat) <- 1e4

mypriors.mod <- PrevMap::control.prior(beta.mean = c(0, 0),
                                       beta.covar = covarsmat,
                                       log.normal.nugget = c(-5,5), # this is tau2
                                       uniform.phi = c(0,10),
                                       log.normal.sigma = c(0,25)
)

mcmcdirections.intercept <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                        n.sim = 1e4+1e3,
                                                        thin = 1, # don't thin
                                                        L.S.lim = c(5,50),
                                                        epsilon.S.lim = c(0.01, 0.1),
                                                        start.nugget = 0.05,
                                                        start.sigma2 = 3,
                                                        start.beta = -5,
                                                        start.phi = 0.5,
                                                        start.S = predict(fit.glm))

mcmcdirections.mod <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                  n.sim = 1e4+1e3,
                                                  thin = 1, # don't thin
                                                  L.S.lim = c(5,50),
                                                  epsilon.S.lim = c(0.01, 0.1),
                                                  start.nugget = 0.05,
                                                  start.sigma2 = 3, 
                                                  start.beta = c(-5, -2),
                                                  start.phi = 0.5,
                                                  start.S = predict(fit.glm))



mod.framework$mypriors <- lapply(1:nrow(mod.framework), function(x) return(mypriors.intercept))
mod.framework$mcmcdirections <- lapply(1:nrow(mod.framework), function(x) return(mcmcdirections.intercept))


# NOTE THIS HACK here, need to have multiple starting betas for the multiple betas in the model
mod.framework$mypriors[[2]] <- mypriors.mod
mod.framework$mcmcdirections[[2]] <- mcmcdirections.mod


#......................
# Make a wrapper for PrevMap
#......................
fit_bayesmap_wrapper <- function(name, 
                                 formula, 
                                 trials = "n", 
                                 coords, 
                                 data, 
                                 mypriors, mcmcdirections, kappa = 1.5){
  
  ret <- PrevMap::binomial.logistic.Bayes(
    formula = as.formula(formula),
    units.m = as.formula(paste("~", trials)),
    coords = coords,
    data = data,
    control.prior = mypriors,
    control.mcmc = mcmcdirections,
    kappa = kappa
  )
  return(ret)
}



# replicate this four times for our four chains
mod.framework <- lapply(1:4, function(x) return(mod.framework)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(name)

####################################################################################
###########                      Diagnostic Runs                           #########
#####################################################################################


setwd("analyses/07-spatial_prediction")
ntry <- nrow(mod.framework)
sjob <- rslurm::slurm_apply(f = fit_bayesmap_wrapper, 
                            params = mod.framework, 
                            jobname = 'Prevmap_Diagnostic_Chains',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))


####################################################################################
###########                          LONG Run                              #########
####################################################################################
# note that PrevMap::spatial.pred.binomial.Bayes access the "formula"
# call within the model object for prediction (with covariates). As a result,
# we can't use a wrapper with purrr as above. 

# Directions LONG RUN                      
mcmcdirections.intercept <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                        n.sim = 1e5+1e3,
                                                        thin = 1, # don't thin
                                                        L.S.lim = c(5,50),
                                                        epsilon.S.lim = c(0.01, 0.1),
                                                        start.nugget = 0.05,
                                                        start.sigma2 = 3,
                                                        start.beta = -5,
                                                        start.phi = 0.5,
                                                        start.S = predict(fit.glm))

mcmcdirections.mod <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                  n.sim = 1e5+1e3,
                                                  thin = 1, # don't thin
                                                  L.S.lim = c(5,50),
                                                  epsilon.S.lim = c(0.01, 0.1),
                                                  start.nugget = 0.05,
                                                  start.sigma2 = 3, 
                                                  start.beta = c(-5, -2),
                                                  start.phi = 0.5,
                                                  start.S = predict(fit.glm))


longrun.prevmapbayes.intercept <- PrevMap::binomial.logistic.Bayes(
  formula = as.formula("plsmdlogit ~ 1"),
  units.m = as.formula("~ n"),
  coords = as.formula("~ longnum + latnum"),
  data = pvclust.weighted.nosf,
  control.prior = mypriors.intercept,
  control.mcmc = mcmcdirections.intercept,
  kappa = 1.5
)

longrun.prevmapbayes.mod <- PrevMap::binomial.logistic.Bayes(
  formula = as.formula("plsmdlogit ~ 1 + precip_mean_cont_scale_clst"),
  units.m = as.formula("~ n"),
  coords = as.formula("~ longnum + latnum"),
  data = pvclust.weighted.nosf,
  control.prior = mypriors.mod,
  control.mcmc = mcmcdirections.mod,
  kappa = 1.5
)

# save out
dir.create("Prevmap_Long_Chains")
saveRDS(object = longrun.prevmapbayes.intercept, file = "Prevmap_Long_Chains/longrun-prevmapbayes-intercept.rds")
saveRDS(object = longrun.prevmapbayes.mod, file = "Prevmap_Long_Chains/longrun-prevmapbayes-mod.rds")


