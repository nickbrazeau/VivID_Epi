library(tidyverse)
library(sf)
library(PrevMap)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")




#----------------------------------------------------------------------------------------------------
# Mapping Functions
#----------------------------------------------------------------------------------------------------
prev_point_est_summarizer <- function(design, maplvl, plsmdmspec, adm1shp){
  
  # rlang
  maplvl <- enquo(maplvl)
  plsmdmspec <- enquo(plsmdmspec)
  
  # clusters are weighted (each individual has same weight in cluster)
  ret <- design %>% 
    dplyr::mutate(count = 1) %>% 
    dplyr::group_by(!!maplvl) %>% 
    dplyr::summarise(n = srvyr::survey_total(count), 
                     plsmdn = srvyr::survey_total(!!plsmdmspec, na.rm = T), 
                     plsmdprev = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)
                     ) %>% 
    dplyr::mutate(logitplsmdprev = logit(plsmdprev, tol = 1e-3))
  
  
  if( quo_name(maplvl) == "adm1name" ){
    ret <- inner_join(ret, adm1shp, by = "adm1name")
   
  } else if( quo_name(maplvl) == "hv001" ){
    clustgeom <- design %>% 
      dplyr::select(c("hv001", "longnum", "latnum", "geometry")) %>% 
      as.data.frame(.) %>% 
      dplyr::filter(!duplicated(.))
    
    ret <- inner_join(ret, clustgeom, by = "hv001")
    
  } else {
    stop("maplvl is not in the options for this function")
  }
  
  # return
  return(ret)
}



mapplotter <- function(data, maplvl, plsmdmspec){
  
  ret <- ggplot() + 
    geom_sf(data = DRCprov) +
    ggtitle(paste(plsmdmspec)) +
    vivid_theme +
    theme(axis.text = element_blank(),
          axis.line = element_blank())
  
  if(maplvl == "adm1name"){
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmdprev)) +
      # scale_fill_distiller("Prevalence", palette = "Spectral") +
      scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
      # scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmdprev[data$plsmdprev != 0], 0.75)) + 
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else if(maplvl == "hv001"){
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmdprev, colour = plsmdprev, size = n), alpha = 0.4) +
      # scale_fill_distiller("Prevalence", palette = "Spectral") +
      scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
      # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmdprev[data$plsmdprev != 0], 0.75)) + 
      scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else {
    stop("maplvl is not in the options for this function")
  }
  
  
  return(ret)
  
}


casemapplotter <- function(data, plsmdmspec){
  # Set some colors ; took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/ ; Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars

  pos <- data %>% 
    dplyr::filter(plsmdprev > 0)
  neg <- data %>% 
    dplyr::filter(plsmdprev == 0)
  
  ret <- ggplot() + 
    geom_sf(data = DRCprov) +
    geom_jitter(data = neg, aes(x=longnum, y=latnum, size = n), shape = 4, show.legend = F, colour = "#377eb8") +
    geom_point(data = pos, aes(x=longnum, y=latnum, colour = plsmdprev, size = n), alpha = 0.4) +
    scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
    scale_size(guide = 'none') +
    ggtitle(paste(plsmdmspec)) +
    coord_sf(datum=NA) + # to get rid of gridlines
    vivid_theme +
    theme(axis.text = element_blank(),
          axis.line = element_blank(), 
          axis.title = element_blank(),
          legend.position = "bottom")
  
  return(ret)
  
}



#----------------------------------------------------------------------------------------------------
# Guassian Surface Plot
#----------------------------------------------------------------------------------------------------
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




prevmaprasterplotter <- function(prevrasters, smoothfct = 5, alpha = 0.8){
  
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
    geom_raster(data = ret.smrstr.m.df, aes(lon, lat, fill = prev), alpha = alpha) +
    scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") 

  return(ret.smrstr.m.plot)
}
