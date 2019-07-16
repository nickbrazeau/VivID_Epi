



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
