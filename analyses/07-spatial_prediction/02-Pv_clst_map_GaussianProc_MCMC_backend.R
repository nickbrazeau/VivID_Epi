#----------------------------------------------------------------------------------------------------
# Purpose of this script is to produce spatial models
# from prevmap for Pv 
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R")
source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) 
library(raster)
library(PrevMap)
set.seed(48, "L'Ecuyer")
tol <- 1e-3
#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
#......................
# Subset to Pv
#......................
pvclust.weighted.nosf <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = srvyr::survey_total(count), 
                   plsmdn = srvyr::survey_total(pv18s, na.rm = T), 
                   plsmdprev = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95))

# need to keep integers, so will round
pvclust.weighted.nosf <- pvclust.weighted.nosf %>% 
  dplyr::mutate(plsmdn = round(plsmdn, 0),
                n = round(n, 0))

#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------
# precipitation already in dataset 
pvclst.covar <- dt %>% 
  dplyr::select("hv001", "precip_mean_cont_scale_clst") %>% 
  dplyr::filter(!duplicated(.))

# bring in cropland
crop <- readRDS("data/derived_data/vividepi_cropland_propmeans.rds") %>% 
  dplyr::select(c("hv001", "cropprop_cont_scale_clst"))
# bring in urban friction
urbanfric <- readRDS("data/derived_data/vividepi_fricurban_clstmeans.rds") %>% 
  dplyr::select(c("hv001", "frctmean_cont_scale_clst"))

pvclst.covar <- dplyr::left_join(x = pvclst.covar, y = crop, by = "hv001") %>% 
  dplyr::left_join(x = ., y = urbanfric, by = "hv001")

# join together
pvclust.weighted.nosf <- dplyr::left_join(pvclust.weighted.nosf, pvclst.covar, by = "hv001")

riskvars = c("precip_mean_cont_scale_clst", 
             "cropprop_cont_scale_clst", "frctmean_cont_scale_clst")



#-------------------------------------------------------------------------
# Make Model Framework
#-------------------------------------------------------------------------
mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                formula = c("plsmdn ~ 1", 
                                            paste0("plsmdn ~ ", paste(riskvars, collapse = "+")))
)

coords <- as.formula(paste0("~ longnum + latnum"))
mod.framework$data <- lapply(1:nrow(mod.framework), function(x) return(pvclust.weighted.nosf))
mod.framework$coords <- lapply(1:nrow(mod.framework), function(x) return(coords))



####################################################################################
###########           PRIORS  & Direction Diagnostic                      ##########
####################################################################################
fit.glm <- glm(cbind(plsmdn, n - plsmdn) ~ 1, 
               data = pvclust.weighted.nosf,
               family = binomial)


mypriors.intercept <- PrevMap::control.prior(beta.mean = 0,
                                             beta.covar = 1,
                                             log.normal.nugget = c(0, 2.5), # this is tau2
                                             uniform.phi = c(0,10),
                                             log.normal.sigma = c(0, 2.5))

# NB covar matrix
covarsmat <- matrix(0, ncol = 4, nrow=4) # three risk factors, 4 betas
diag(covarsmat) <- 1 # identity matrix

mypriors.mod <- PrevMap::control.prior(beta.mean = c(0, 0, 0, 0),
                                       beta.covar = covarsmat,
                                       log.normal.nugget = c(0, 2.5), # this is tau2
                                       uniform.phi = c(0,10),
                                       log.normal.sigma = c(0, 2.5))

mcmcdirections.intercept <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                        n.sim = 1e4+1e3,
                                                        thin = 10, 
                                                        L.S.lim = c(5,50),
                                                        epsilon.S.lim = c(0.01, 0.1),
                                                        start.nugget = 1,
                                                        start.sigma2 = 0.2,
                                                        start.beta = -5,
                                                        start.phi = 0.5,
                                                        start.S = predict(fit.glm))

mcmcdirections.mod <- PrevMap::control.mcmc.Bayes(burnin = 1e3, 
                                                  n.sim = 1e4+1e3,
                                                  thin = 10, 
                                                  L.S.lim = c(5,50),
                                                  epsilon.S.lim = c(0.01, 0.1),
                                                  start.nugget = 1,
                                                  start.sigma2 = 0.2, 
                                                  start.beta = c(-4, -0.2, 0, -0.02),
                                                  start.phi = 0.5,
                                                  start.S = predict(fit.glm))


mod.framework$mypriors <- lapply(1:nrow(mod.framework), 
                                 function(x) return(mypriors.intercept))
mod.framework$mcmcdirections <- lapply(1:nrow(mod.framework), 
                                       function(x) return(mcmcdirections.intercept))


# NOTE THIS HACK here, need to have multiple starting betas for the multiple betas in the model
mod.framework$mypriors[[ which(stringr::str_detect(mod.framework$formula, "\\+")) ]] <- mypriors.mod
mod.framework$mcmcdirections[[ which(stringr::str_detect(mod.framework$formula, "\\+")) ]] <- mcmcdirections.mod


#......................
# Make a wrapper for PrevMap
#......................
fit_bayesmap_wrapper <- function(name, 
                                 formula, 
                                 trials = "n", 
                                 coords, 
                                 data, 
                                 mypriors, mcmcdirections, 
                                 kappa = 0.75){
  
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
                            jobname = 'prevmap_diagnostic_chains',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 64000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5-00:00:00"))

####################################################################################
###########                          LONG Run                              #########
####################################################################################
# note that PrevMap::spatial.pred.binomial.Bayes access the "formula"
# call within the model object for prediction (with covariates). As a result,
# we can't use a wrapper with purrr as above. 

# Directions LONG RUN                      
mcmcdirections.intercept <- PrevMap::control.mcmc.Bayes(burnin = 1e4, 
                                                        n.sim = 1e5 + 1e4,
                                                        thin = 100, 
                                                        L.S.lim = c(5,50),
                                                        epsilon.S.lim = c(0.01, 0.1),
                                                        start.nugget = 1,
                                                        start.sigma2 = 0.2,
                                                        start.beta = -5,
                                                        start.phi = 0.5,
                                                        start.S = predict(fit.glm))

mcmcdirections.mod <- PrevMap::control.mcmc.Bayes(burnin = 1e4, 
                                                  n.sim = 1e5 + 1e4,
                                                  thin = 100, # don't thin
                                                  L.S.lim = c(5,50),
                                                  epsilon.S.lim = c(0.01, 0.1),
                                                  start.nugget = 1,
                                                  start.sigma2 = 0.2, 
                                                  start.beta = c(-4, -0.2, 0, -0.02),
                                                  start.phi = 0.5,
                                                  start.S = predict(fit.glm))



longrun.prevmapbayes.intercept <- PrevMap::binomial.logistic.Bayes(
  formula = as.formula("plsmdn ~ 1"),
  units.m = as.formula("~ n"),
  coords = as.formula("~ longnum + latnum"),
  data = pvclust.weighted.nosf,
  control.prior = mypriors.intercept,
  control.mcmc = mcmcdirections.intercept,
  kappa = 0.75
)


longrun.prevmapbayes.mod <- PrevMap::binomial.logistic.Bayes(
  formula = as.formula("plsmdn ~ 1 + precip_mean_cont_scale_clst + 
                       cropprop_cont_scale_clst + frctmean_cont_scale_clst"),
  units.m = as.formula("~ n"),
  coords = as.formula("~ longnum + latnum"),
  data = pvclust.weighted.nosf,
  control.prior = mypriors.mod,
  control.mcmc = mcmcdirections.mod,
  kappa = 0.75
)



# save out
dir.create("prevmap_long_chains")
saveRDS(object = longrun.prevmapbayes.intercept, file = "prevmap_long_chains/longrun-prevmapbayes-intercept.rds")
saveRDS(object = longrun.prevmapbayes.mod, file = "prevmap_long_chains/longrun-prevmapbayes-mod.rds")


