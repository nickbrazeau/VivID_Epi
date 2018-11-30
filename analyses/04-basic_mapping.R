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




#----------------------------------------------------------------------------------------------------
# Explore here the different prevalences/regions of Plasmodium species
#----------------------------------------------------------------------------------------------------
# issue here is that we are not using the weights. easy to incorporate...but is it worthwhile
# at the very least, need to include the HIV weights

prev_summarizer <- function(x = dtsub, provlvl, plsmdmspec, sfobj){
  if(is.character(x)){
    stop("this is the issue")
  }
  # catch error
  # if( !( any(colnames(x) %in% c("hv001")) & any(colnames(x) %in% c("hv005")) ) ){
  #   stop("Must include cluster and cluster weight. need to check if you are using pr or hiv sample weights")
  # }
  # rlang
  provlvl <- enquo(provlvl)
  plsmdmspec <- enquo(plsmdmspec)
  
  # clusters are weighted (each individual has same weight in cluster)
   ret <- x %>% 
    srvyr::as_survey_design(ids = hv001, weights = hv005) %>% 
    dplyr::group_by(!!provlvl) %>% 
    dplyr::summarise(plsmd = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)
    ) %>% 
     left_join(., sfobj) # attach spatial data, let R figure out the common var
   # return
   return(ret)
}

# set mapping parameters
map_params <- expand.grid(
  plsmdmspec = c("pfldh", "pv18s", "po18s"),
  provlvl = c("adm1name", "adm1name"), 
  stringsAsFactors = FALSE
)


map_params$sfobj <- ifelse(map_params$provlvl == "adm1name", "DRCprov", 
                           ifelse(map_params$provlvl == "hv001", "ge", 
                                  NA))


purrr::pmap(.l = map_params, .f = prev_summarizer)

prev_summarizer(x = dtsub, provlvl = hv001, plsmdmspec = pfldh, sfobj = ge)
prev_summarizer(x = dtsub, provlvl = adm1name, plsmdmspec = pfldh, sfobj = DRCprov)

p

x %>% 
  ggplot(data = .) + 
  geom_sf(aes(fill = plsmd)) 



# use purrr ot attach maps and plot

drugrx_wsaf_df <- drugrx_wsaf %>%
  tidyr::gather(., key = "Sample_ID", value = "wsaf", 4:ncol(.)) %>%
  dplyr::left_join(x=p$data_processed$data_sample, y=., by = "Sample_ID") %>%
  dplyr::group_by(Country, name) %>%
  tidyr::nest() %>%
  dplyr::mutate(datafilt = purrr::map(data, smplfilter, cov = 0.5))


drugrx_wsaf_df <- drugrx_wsaf_df %>%
  mutate(plotObj = purrr::pmap(list(df = drugrx_wsaf_df$datafilt,
                                    site = drugrx_wsaf_df$Country,
                                    drugrx_exp = lapply(1:length(drugrx_wsaf_df$datafilt), function(x)return(drugrx_exp)),
                                    loci = drugrx_wsaf_df$name),
                               plotwsaf)
  )










# Set some colors ; took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/ ; Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars
prevscale <- rev(heat.colors(101))
cool <- rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm <- rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols <- c(rev(cool), rev(warm))
mypalettediff <- colorRampPalette(cols)(101)

 
 
    
    dhsprevbyprov07 <- dhs %>% 
    dplyr::select(provname, cluster, pf, year, long, lat) %>% 
    dplyr::filter(year == "2007") %>% 
    dplyr::filter(!is.na(cluster)) %>% 
    dplyr::group_by(provname) %>% 
    dplyr::summarise(meanpf=mean(pf, na.rm=T), meanlat=mean(lat), meanlong=mean(long))
  
  FALSE %in% c(DRCprov@data$REGNAME %in% dhs$provname) # are prov names same
  
  DRCprov07 <- sp::merge(DRCprov, dhsprevbyprov07, by.x="REGNAME", by.y="provname")
  DRCprov07plot <- sp::spplot(DRCprov07, "meanpf", main="DHS 2007 Prevalence by Province", col.regions = prevscale, cuts = 100)
}






#----------------------------------------------------------------------------------------------------
# Purpose of this script is to explore basic maps of Plasmodium infections in the CD2013 data
#----------------------------------------------------------------------------------------------------
# add in GPS points
cd2013gps <- tibble(cluster = cd2013dhsge$DHSCLUST, lat = cd2013dhsge$LATNUM, long = cd2013dhsge$LONGNUM)
# Clusters with Lat and Long of 0,0 were not able to be identified and should have coordinates set to NA
cd2013gps <- cd2013gps %>%
  dplyr::mutate(long = ifelse(long == 0 & lat == 0, NA, long)) %>%
  dplyr::mutate(lat = ifelse(long == 0 & lat == 0, NA, lat))
