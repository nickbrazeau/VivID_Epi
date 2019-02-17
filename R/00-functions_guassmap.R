library(tidyverse)
library(sf)
library(PrevMap)
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions_basic.R")

#----------------------------------------------------------------------------------------------------
# Guassian Surface Plot
#----------------------------------------------------------------------------------------------------

load("~/Documents/GitHub/VivID_Epi/data/vividmaps_small.rda")
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum")]


guass_map_clstr_summarizer <- function(data, plsmdmspec){
  
  clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum")]
  
  plsmdmspec <- enquo(plsmdmspec)
  # clusters are weighted (each individual has same weight in cluster)
  ret <- data %>% 
    dplyr::mutate(count = 1) %>% 
    srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
    dplyr::group_by(hv001) %>% 
    dplyr::summarise(n = srvyr::survey_total(count), 
                     plsmdn = srvyr::survey_total(!!plsmdmspec, na.rm = T), 
                     plsmdprev = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)
                     ) %>% 
    dplyr::inner_join(., clustgeom, by = "hv001")
  ret <- ret %>% 
    dplyr::mutate(logitplsmdprev = logit(plsmdprev)) # weird how this is being held
  
  # return
  return(ret)
}



fit_pred_spMLE <- function(outcome, covar, 
                     long_var = "longnum", lat_var = "latnum",
                     data, 
                     grid.pred,
                     kappa = 0.5,
                     start.cov.pars = c(1,1),
                     scale.predictions = "prevalence",
                     pred.reps = 1e2, SE = T){
  eq <- as.formula(paste0(outcome, "~", covar))
  coords <- as.formula(paste0("~", long_var, "+", lat_var))
  ret.fit <- PrevMap::linear.model.MLE(formula=eq, coords=coords, 
                          data=data, start.cov.pars=start.cov.pars, 
                          kappa=kappa)
  
  ret.pred <- PrevMap::spatial.pred.linear.MLE(ret.fit, 
                          grid.pred = grid.pred, 
                          scale.predictions=scale.predictions, 
                          n.sim.prev=pred.reps, standard.errors=SE)
  
  return(list(
    fit = ret.fit,
    pred = ret.pred
  ))
  
}




prevmaprasterplotter <- function(prevrasters, smoothfct = 5){
  
  ret.rstr <- raster::rasterFromXYZ(cbind(prevrasters$grid.pred[,1],
                                          prevrasters$grid.pred[,2],
                                          prevrasters$prevalence$predictions),
                                      crs="+proj=longlat +datum=WGS84")
  
  ret.smrstr <- focal(ret.rstr, w=matrix(1,
                                               nrow=smoothfct,
                                               ncol=smoothfct), mean)
  
  ret.smrstr.m  <-  raster::rasterToPoints(ret.smrstr)
  ret.smrstr.m.df <-  data.frame(lon = ret.smrstr.m[,1], 
                                 lat = ret.smrstr.m[,2], 
                                 prev = ret.smrstr.m[,3])
  ret.smrstr.m.plot <- ggplot() + 
    geom_raster(data = ret.smrstr.m.df, aes(lon, lat, fill = prev), alpha = 0.8) +
    scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") 

  return(ret.smrstr.m.plot)
}
