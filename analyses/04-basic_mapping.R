#----------------------------------------------------------------------------------------------------
# Purpose of this script is to explore basic maps of Plasmodium infections in the CD2013 data
#----------------------------------------------------------------------------------------------------
source("analyses/00-functions.R") 
library(tidyverse)
library(sf)
library(srvyr) #wrap the survey package in dplyr syntax
#......................
# Import Data
#......................
load("data/vividepi_raw.rda")
dt <- merge_pr_plsmdm_gemtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge)



#......................
# Datawrangle
#......................
dtsub <- dt %>% 
  dplyr::select(c("hivrecode_barcode", "hv001", "hv002", "hv005", "hiv05", "pfldh", "po18s", "pv18s", "adm1dhs", "adm1name", "latnum", "longnum")) %>% 
  dplyr::mutate(hv005 = hv005/1e5,
                hiv05 = hiv05/1e5) # per the DHS website https://dhsprogram.com/data/Using-DataSets-for-Analysis.cfm#CP_JUMP_14042 --- like other variables in DHS datasets, decimal points are not included in the weight variable. Analysts need to divide the sampling weight they are using by 1,000,000
#   tidyr::gather(., key="plsmdmspec", value = "outcome", 6:8) %>%
#   dplyr::group_by(plsmdmspec) %>% 
#   tidyr::nest() 
# 
# 
# dtsub_summ <- dplyr::bind_cols(dplyr::bind_rows(dtsub, dtsub), maplvl = c(rep("adm1name", nrow(dtsub)), rep("hv001", nrow(dtsub)))) %>% 
#   dplyr::mutate(sfobj = ifelse(maplvl == "adm1name", "DRCprov", 
#                                  ifelse(maplvl == "hv001", "ge", 
#                                         NA))
#   )

#----------------------------------------------------------------------------------------------------
# Explore here the different prevalences/regions of Plasmodium species basic maps
#----------------------------------------------------------------------------------------------------
# issue here is that we are not using the weights. easy to incorporate...but is it worthwhile
# at the very least, need to include the HIV weights

prev_summarizer <- function(data, maplvl, plsmdmspec, sfobj){
  
  
  # catch error
  # if( !( any(colnames(x) %in% c("hv001")) & any(colnames(x) %in% c("hv005")) ) ){
  #   stop("Must include cluster and cluster weight. need to check if you are using pr or hiv sample weights")
  # }
  # rlang
  maplvl <- enquo(maplvl)
  plsmdmspec <- enquo(plsmdmspec)
  
  # clusters are weighted (each individual has same weight in cluster)
  ret <- data %>% 
    srvyr::as_survey_design(ids = hv001, weights = hv005) %>% 
    dplyr::group_by(!!maplvl) %>% 
    dplyr::summarise(plsmdn = srvyr::survey_total(hv001), 
                     plsmd = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)) %>% 
    dplyr::left_join(., sfobj) # attach spatial data, let R figure out the common var
  # return
  return(ret)
}


# set mapping parameters
# map_params <- expand.grid(
#   plsmdmspec = c("pfldh", "pv18s", "po18s"),
#   provlvl = c("adm1name", "hv001"), 
#   stringsAsFactors = FALSE
# )
# 
# map_params$sfobj <- ifelse(map_params$provlvl == "adm1name", "DRCprov", 
#                            ifelse(map_params$provlvl == "hv001", "ge", 
#                                   NA))
# NOTE PURRR is having an issue because of the "~"


pfldhprov <- prev_summarizer(data = dtsub, maplvl = adm1name, plsmdmspec = pfldh, sfobj = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "adm1name")
pv18sprov <- prev_summarizer(data = dtsub, maplvl = adm1name, plsmdmspec = pv18s, sfobj = DRCprov)  %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "adm1name")
po18sprov <- prev_summarizer(data = dtsub, maplvl = adm1name, plsmdmspec = po18s, sfobj = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "adm1name")

pfldhclust <- prev_summarizer(data = dtsub, maplvl = hv001, plsmdmspec = pfldh, sfobj = ge) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001")
pv18sclust <- prev_summarizer(data = dtsub, maplvl = hv001, plsmdmspec = pv18s, sfobj = ge) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001")
po18sclust <- prev_summarizer(data = dtsub, maplvl = hv001, plsmdmspec = po18s, sfobj = ge) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001")


# this awful hack becuase of this issue https://github.com/tidyverse/dplyr/issues/3483
# we are going down the rabbit hole just to try and make this stupid survey and purr package work. fine for now but return
mp <- dplyr::bind_rows(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust) %>% 
  dplyr::group_by(plsmdmspec, maplvl) %>% 
  tidyr::nest()

mp$data <- sapply(list(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust), function(x) return(x))


#......................
# Plot Maps
#......................
mapplotter <- function(data, maplvl, plsmdmspec){
  
  
  ret <- ggplot() + 
    geom_sf(data = DRCprov) +
    ggtitle(paste(plsmdmspec)) +
    theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 13), 
          legend.position = "bottom",
          legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
          legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
          axis.text = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.grid = element_blank(),
          panel.border = element_blank())
  
  if(maplvl == "adm1name"){
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmd)) +
      scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else if(maplvl == "hv001"){
    
    ret <- ret + geom_sf(data = data, aes(fill = plsmd, colour = plsmd, size = plsmdn), alpha = 0.8) +
      scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
      scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
      coord_sf(datum=NA)  # to get rid of gridlines
    
  } else {
    stop("maplvl is not in the options for this function")
  }
  
  
  return(ret)
  
}

prevmaps <- pmap(mp, mapplotter)

jpeg(file = "figures/04-prevmaps.jpg", width = 11, height = 8, units = "in", res=300)
gridExtra::grid.arrange(prevmaps[[1]], 
                        prevmaps[[2]],
                        prevmaps[[3]],
                        prevmaps[[4]],
                        prevmaps[[5]],
                        prevmaps[[6]],
                        nrow=2, top=grid::textGrob("Prevalence by Species in CD2013 DHS", gp=grid::gpar(fontsize=15, fontfamily = "Arial", fontface = "bold"))) 
graphics.off()

#----------------------------------------------------------------------------------------------------
# Explore here the different clusters of Plasmodium species z scores

clusters <- mp %>% 
  dplyr::filter(maplvl == "hv001")

calczscoreprev <- function(data){
  data <- data %>% 
    dplyr::mutate(meanplsmd = mean(plsmd, na.rm=T),
                  sdplsmd = sd(plsmd, na.rm=T),
                  zscore = (plsmd - meanplsmd)/sdplsmd
                  )
  return(data)
}


clusters$ret <- map(clusters$data,
                     calczscoreprev)

#......................
# Plot Maps
#......................
zscoremapplotter <- function(data, plsmdmspec){
  # Set some colors ; took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/ ; Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars
  
  ret <- data %>% 
    ggplot() + 
    geom_sf(data = DRCprov) +
    geom_sf(data = data, aes(colour = zscore, size = plsmdn), alpha = 0.6) +
    scale_color_gradient2("Prevalence Z-Scores", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = median(data$plsmd[data$plsmd != 0])) + 
    scale_size(guide = 'none') +
    ggtitle(paste(plsmdmspec)) +
    coord_sf(datum=NA) + # to get rid of gridlines
    theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 18), 
          legend.position = "bottom",
          legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
          legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
          axis.text = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.grid = element_blank(),
          panel.border = element_blank())
  
  return(ret)
  
}

zscoreprevmaps <- pmap(list(data = clusters$ret, plsmdmspec = clusters$plsmdmspec), zscoremapplotter)

prevhist <- pmap(list(data = clusters$ret, plsmdmspec = clusters$plsmdmspec), function(data, plsmdmspec){
  data %>% 
    dplyr::filter(plsmd != 0) %>% 
  ggplot(.) +
    geom_histogram(mapping = aes(x=plsmd)) + 
    ggtitle(paste(plsmdmspec)) + ylab("Count") +
    theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 18), 
          legend.position = "bottom",
          legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
          legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 60, size = 10),
          axis.text.y = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, size = 10),
          axis.ticks = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          panel.grid = element_blank(),
          panel.border = element_blank())
})



jpeg(file = "figures/04-zscoremaps.jpg", width = 11, height = 8, units="in", res=300)
gridExtra::grid.arrange(prevhist[[1]],
                        prevhist[[2]],
                        prevhist[[3]],
                        zscoreprevmaps[[1]], 
                        zscoreprevmaps[[2]],
                        zscoreprevmaps[[3]],
                        nrow=2, 
                        top=grid::textGrob("Histograms and Standardized Prevalence by Species in CD2013 DHS", gp=grid::gpar(fontsize=15, fontfamily = "Arial", fontface = "bold"))) 
graphics.off()

