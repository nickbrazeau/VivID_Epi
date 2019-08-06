#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create a spatial prediction model
# from PrevMap that is based on a Gaussain Process Model
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) 
library(PrevMap)
set.seed(48)
tol <- 1e-3
#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
dtsrvy <- makecd2013survey(survey = dt)
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
ge <- sf::st_as_sf( readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds") )


#......................
# Subset to Pv
#......................
pvclust.weighted <- mp %>% 
  dplyr::filter(plsmdmspec == "pv18s" & maplvl == "hv001") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pvclust.weighted <- sf::st_as_sf(pvclust.weighted)
sf::st_crs(pvclust.weighted) <-  sf::st_crs(ge)

# need to keep integers, so will round
pvclust.weighted <- pvclust.weighted %>% 
  dplyr::mutate_if(is.numeric, round, 0)

#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------

pvclst.covar <- dtsrvy %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(meanprov_precip_lag_cont_scale_clst = srvyr::survey_mean(precip_ann_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95),
                   meanprov_alt_dem_cont_scale_clst = srvyr::survey_mean(alt_dem_cont_scale_clst, na.rm = T, vartype = c("se", "ci"), level = 0.95)
  )

pvclust.weighted <- dplyr::left_join(pvclust.weighted, pvclst.covar, by = "hv001")
pvclust.weighted.nosf <- pvclust.weighted
sf::st_geometry(pvclust.weighted.nosf) <- NULL


#----------------------------------------------------------------------------------------------------
# PrevMap & The Gaussian Process Model 
#----------------------------------------------------------------------------------------------------
# Following Giorgi & Diggle 2017 -- https://www.jstatsoft.org/article/view/v078i08
#......................
# EDA - Variogram
#......................
eda <- pvclust.weighted.nosf %>% 
  dplyr::select(c("hv001", "plsmdn", "n", "plsmdprev", "longnum", "latnum"))

# transform count of "successes" to logit space
eda$plsmdnlogit <- log( (eda$plsmdn + 0.5)/(eda$n - eda$plsmdn + 0.5) ) # 0.5 as tolerance for 0s

# look at kappa smooth parameter for matern covariance
profile.kappa <- PrevMap::shape.matern(plsmdnlogit ~ 1,
                                       coords = ~ longnum + latnum,
                                       data = eda,
                                       set.kappa = seq(0.1, 10, length = 16),
                                       start.par = c(0.2, 0.05), # starting values for the scale parameter phi and the relative variance of the nugget effect nu2
                                       coverage = NULL # CI coverage
                                      )


# fit a voariogram 
coords <- eda[,c("longnum", "latnum")]
# note, this defaults to euclidean distance (from coordinates of long,lat)
vari <- geoR::variog(coords = coords, data = eda$plsmdnlogit, max.dist = 5)
vari.fit <- geoR::variofit(vari, 
                           ini.cov.pars = c(3, 1),
                           cov.model = "matern",
                           fix.kappa = F)


plot(vari, main = substitute(paste("Semivariogram for ", italic("P. vivax"))))
# really not an impressive amount of spatial autocorrelation

#......................
# Make Model Framework
#......................
clst.covar.names <- c("precip_ann_cont_scale_clst", "meanprov_alt_dem_cont_scale_prov")
mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                formula = c("plsmdn ~ 1", 
                                            paste0("plsmdn ~ ", paste(clst.covar.names, collapse = " + "))),
                                family = "binomial"
                                )

coords <- as.formula(paste0("~ longnum + latnum"))
trials <- pvclust.weighted$n

#......................
# PRIORS
#......................
mypriors <- PrevMap::control.prior()
mcmcdirections <- PrevMap::control.mcmc.Bayes(burnin = 1e1, 
                                              n.sim = 1e2,
                                              thin = 0
                                              )

#......................
# Make a wrapper for PrevMap
#......................
wrap_S.binomial.logistic.Bayes <- function(...){
  PrevMap::binomial.logistic.Bayes(
                                   formula = formula,
                                   units.m = trials, 
                                   coords = coords, 
                                   data = data, 
                                   
  
  
)

(formula=eq, coords=coords, 
                                     data=data, start.cov.pars=start.cov.pars, 
                                     kappa=kappa)




  
  
  
  
  
  
  
  
  
  
  set.seed(1234)
  data(data_sim)
  # Select a subset of data_sim with 50 observations
  n.subset <- 50
  data_subset <- data_sim[sample(1:nrow(data_sim),n.subset),]
  # Set the MCMC control parameters
  control.mcmc <- control.mcmc.Bayes(n.sim=10,burnin=0,thin=1,
                                     h.theta1=0.05,h.theta2=0.05,
                                     L.S.lim=c(1,50),epsilon.S.lim=c(0.01,0.02),
                                     start.beta=0,start.sigma2=1,start.phi=0.15,
                                     start.nugget=NULL,
                                     start.S=rep(0,n.subset))
  
  cp <- control.prior(beta.mean=0,beta.covar=1,
                      log.normal.phi=c(log(0.15),0.05),
                      log.normal.sigma2=c(log(1),0.1))
  
  fit.Bayes <- binomial.logistic.Bayes(formula=y~1,coords=~x1+x2,units.m=~units.m,
                                       data=data_subset,control.prior=cp,
                                       control.mcmc=control.mcmc,kappa=2)
  summary(fit.Bayes)
  par(mfrow=c(2,4))
  autocor.plot(fit.Bayes,param="S",component.S="all")
  autocor.plot(fit.Bayes,param="beta",component.beta=1)
  autocor.plot(fit.Bayes,param="sigma2")
  autocor.plot(fit.Bayes,param="phi")
  trace.plot(fit.Bayes,param="S",component.S=30)
  trace.plot(fit.Bayes,param="beta",component.beta=1)
  trace.plot(fit.Bayes,param="sigma2")
  trace.plot(fit.Bayes,param="phi")


























#.............................
# get prev rasters for individual levels
#..............................
# polybb <- osmdata::getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long","lat")

pr$prevrasters <- map(pr$data, 
                      fit_pred_spMLE, outcome = "logitplsmdprev", covar = "1", 
                      long_var = "longnum", lat_var = "latnum",
                      grid.pred = grid.pred, kappa = 0.5, 
                      pred.reps = 10)

pr$prevrasterspred <- purrr::map(pr$prevrasters, "pred")



#.............................
# plot prev rasters
#..............................
prevmaprasterplots <- lapply(pr$prevrasterspred,
                             prevmaprasterplotter, smoothfct = rep(7,3), alpha = 0.5)


#----------------------------------------------------------------------------------------------------
# Save Objects & Write out
#----------------------------------------------------------------------------------------------------

saveRDS(pr, file = "data/derived_data/basic_prevmap_mapping_data.rds")
saveRDS(prevmaprasterplots, file = "results/prevmap_raster_plots.rds")


