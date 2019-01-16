
library(tidyverse)
library(sf)
library(srvyr) # wrap the survey package in dplyr syntax

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
# Tidy and DataWrangling functions
#..................................
merge_pr_plsmdm_gemtdt <- function(pr = arpr, plsmdm = panplasmpcrres, ge = ge){
  
  ret <- dplyr::left_join(plsmdm, pr, by="hivrecode_barcode") %>%
    dplyr::left_join(., ge, by = "hv001")
  
  return(ret)
}



#----------------------------------------------------------------------------------------------------
# Mapping Functions
#----------------------------------------------------------------------------------------------------
prev_point_est_summarizer <- function(data, maplvl, plsmdmspec, sfobj){
  
  
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
    srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
    dplyr::group_by(!!maplvl) %>% 
    dplyr::summarise(plsmdn = srvyr::survey_total(count), 
                     plsmd = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)) %>% 
    dplyr::left_join(., sfobj) # attach spatial data, let R figure out the common var
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
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmd)) +
      scale_fill_distiller("Prevalence", palette = "Spectral") +
     # scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
     # scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else if(maplvl == "hv001"){
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmd, colour = plsmd, size = plsmdn), alpha = 0.8) +
      scale_fill_distiller("Prevalence", palette = "Spectral") +
    #  scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
    # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
      scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else {
    stop("maplvl is not in the options for this function")
  }
  
  
  return(ret)
  
}


mapplotter_clust_terrain <- function(data, plsmdmspec){
  
   ret <-  ggmap::ggmap(drc_stamen_back_terrain) + 
           geom_point(data=data, aes(x=longnum, y = latnum, fill = plsmd, colour = plsmd, size = plsmdn), 
                      alpha = 0.4) +
           scale_fill_distiller("Prevalence", palette = "Spectral") +
          # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
           # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
           scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
           coord_sf(datum=NA) +
           ggtitle(paste(plsmdmspec)) +
           vivid_theme +
           theme(axis.text = element_blank(),
                 axis.line = element_blank(), 
                 legend.position = "bottom")
   
    return(ret)
}
    

