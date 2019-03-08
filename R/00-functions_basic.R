library(tidyverse)
library(sf)
library(srvyr) # wrap the survey package in dplyr syntax


#----------------------------------------------------------------------------------------------------
# Basic epi
#----------------------------------------------------------------------------------------------------

logit <- function(x, tol=1e-4){ 
    return( log(((x+tol)/(1-x+tol))) )
}

#----------------------------------------------------------------------------------------------------
# ViVID Epi Theme
#----------------------------------------------------------------------------------------------------

vivid_theme <- theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 14), 
                     axis.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 12), 
                     axis.text = element_text(famil = "Arial", hjust = 0.5, size = 11), 
                     legend.position = "bottom",
                     legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
                     legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"),
                     panel.grid = element_blank(),
                     panel.border = element_blank())



#----------------------------------------------------------------------------------------------------
# Mapping Functions
#----------------------------------------------------------------------------------------------------
prev_point_est_summarizer <- function(data, maplvl, plsmdmspec){
  
  
  # catch error
  # if( !( any(colnames(x) %in% c("hv001")) & any(colnames(x) %in% c("hv005")) ) ){
  #   stop("Must include cluster and cluster weight. need to check if you are using pr or hiv sample weights")
  # }
  # rlang
  maplvl <- enquo(maplvl)
  plsmdmspec <- enquo(plsmdmspec)
  
  # clusters are weighted (each individual has same weight in cluster)
  ret <- data %>% 
    dplyr::mutate(count = 1) %>% 
    srvyr::as_survey_design(ids = hv001, weights = hiv05_wi) %>% 
    dplyr::group_by(!!maplvl) %>% 
    dplyr::summarise(n = srvyr::survey_total(count), 
                     plsmdn = srvyr::survey_total(!!plsmdmspec, na.rm = T), 
                     plsmdprev = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95))
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
    data <- inner_join(data, DRCprov, by = "adm1name")
    ret <- ret + geom_sf(data = data, aes(fill = plsmdprev)) +
      # scale_fill_distiller("Prevalence", palette = "Spectral") +
       scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
      # scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmdprev[data$plsmdprev != 0], 0.75)) + 
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else if(maplvl == "hv001"){
    clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "geometry")]
    data <- inner_join(data, clustgeom, by = "hv001")

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


mapplotter_clust_terrain <- function(data, plsmdmspec){
  # need to load data/vividmaps_large.rda first and recode.rda
  clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "longnum", "latnum", "geometry")]
  data <- inner_join(data, clustgeom, by = "hv001")
  
  ret <- ggplot() + 
    prettybasemap_terraincolors + 
    geom_point(data = data, aes(x=longnum, y = latnum, colour = plsmdprev, size = n), 
               alpha = 0.4) + # using point here to not overwrite sf fill
    # scale_fill_distiller("Prevalence", palette = "Spectral") +
     scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
    # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmdprev[data$plsmdprev != 0], 0.75)) + 
    scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
    coord_sf(datum=NA) +
    ggtitle(paste(plsmdmspec)) +
    vivid_theme +
    theme(axis.text = element_blank(),
          axis.line = element_blank(), 
          legend.position = "bottom")
  
  return(ret)
}
