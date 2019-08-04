library(tidyverse)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")


# Set some colors
# took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/
# Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars
# color for prev
prevscale <- rev(heat.colors(101))

# diff colors
cool <- rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm <- rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols <- c(rev(cool), rev(warm))
mypalettediff <- colorRampPalette(cols)(101)




#----------------------------------------------------------------------------------------------------
# Read Raster Temp and Precip
#----------------------------------------------------------------------------------------------------
readRasterBB <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = sf::as_Spatial(bb))
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m
  
  return(ret)
  
}

extractDHSpts_from_raster <- function(hv001, raster, 
                                      geometry, urban_rura){
  # note raster is expected to have units of m
  
  geometry <- sf::as_Spatial(geometry)
  buffer = ifelse(urban_rura == "R", 10, 2)
  ret <- raster::extract(x = raster[[1]],
                         y =  geometry[[1]], 
                         buffer = buffer,
                         fun = mean,
                         sp = F) 
  return(ret)
}



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
  
  # get sf class back
  ret <- sf::st_as_sf(ret)
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


casemap_prev_plotter <- function(data, plsmdmspec){
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



casemap_n_plotter <- function(data, plsmdmspec){
  # Set some colors ; took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/ ; Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars
  
  pos <- data %>% 
    dplyr::filter(plsmdn > 0)
  neg <- data %>% 
    dplyr::filter(plsmdn == 0)
  
  ret <- ggplot() + 
    geom_sf(data = DRCprov) +
    geom_jitter(data = neg, aes(x=longnum, y=latnum, size = n), shape = 4, show.legend = F, colour = "#377eb8") +
    geom_point(data = pos, aes(x=longnum, y=latnum, colour = plsmdn, size = n), alpha = 0.4) +
    scale_color_gradient2("Absolute \n Count", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
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


