#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle weather data 
# that is around the time of our study period for the CD2013
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(sp)
library(sf)

#..............................
# Housekeeping
#..............................
readRasterBB.precip <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = bb)
  
  vals <- raster::values(ret) 
  vals <- ifelse(vals == -9999, NA, vals) # improper values
  raster::values(ret) <- vals
  
  
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m
  
  return(ret)
  
}

readRasterBB.temp <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  
  vals <- raster::values(ret) 
  vals <- ifelse(vals <= 7500, NA, vals) # improper values
  vals <- (vals * 0.02) - 273.15
  raster::values(ret) <- vals
  
  ret <- raster::crop(x = ret, y = bb)
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m
  
  return(ret)
  
}


# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"


#......................................................................................................
# Precipitation (CHRIPS) and Temperature (MODIS/LAADS) Read In Data
#......................................................................................................

precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T, 
                     pattern = ".tif")
precipfrst <- lapply(precip, readRasterBB.precip, bb = caf)

precipdf <- tibble::tibble(names = basename(precip)) %>% 
  dplyr::mutate(names = gsub("chirps-v2.0.|.tif", "", names),
                year = stringr::str_split_fixed(names, "\\.", n=2)[,1] ,
                mnth =  stringr::str_split_fixed(names, "\\.", n=2)[,2] ,
                hvdate_dtdmy = lubridate::dmy(paste(1, mnth, year, sep = "/")),
                year = lubridate::year(hvdate_dtdmy),
                mnth = lubridate::month(hvdate_dtdmy),
                hvyrmnth_dtmnth = factor(paste(year, mnth, sep = "-")),
                precipraster = precipfrst) %>% 
  dplyr::select(c("hvyrmnth_dtmnth", "precipraster"))



# NOTE, reading in masked temperature files
tempfiles <- list.files(path = "data/raw_data/weather_data/LAADS_NASA/", full.names = T, 
                       pattern = "LST_Day_CMG.tif")
tempfrst <- lapply(tempfiles, readRasterBB.temp, bb = caf)

tempdf <- tibble::tibble(namestemp = basename(tempfiles)) %>% 
  dplyr::mutate(namestemp = stringr::str_extract(string = namestemp, pattern = "A[0-9][0-9][0-9][0-9][0-9][0-9][0-9]"),
                namestemp = gsub("A", "", namestemp),
                year = substr(namestemp, 1, 4),
                day = substr(namestemp, 5, 7),
                day = as.numeric(day),
                hvdate_dtdy = as.Date(paste0(year, "-", day), format = "%Y-%j", origin = "01-01-2013"),
                hvdate_dtdy = lubridate::ymd(hvdate_dtdy),
                year =  lubridate::year(hvdate_dtdy),
                month = lubridate::month(hvdate_dtdy),
                hvyrmnth_dtmnth = factor(paste(year, month, sep = "-")),
                tempraster = tempfrst
  ) %>% 
  dplyr::select(c("hvyrmnth_dtmnth", "tempraster"))




#......................................................................................................
# Precipitation and Temperature Considered by Month 
#......................................................................................................
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# drop observations with missing geospatial data 
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 

wthrnd.mnth <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth", "geometry", "urban_rura"))


wthrnd.mnth <- wthrnd.mnth %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd.mnth <- wthrnd.mnth[!duplicated(wthrnd.mnth$hv001),]

wthrnd.mnth <- wthrnd.mnth %>% 
  dplyr::left_join(., tempdf) %>% 
  dplyr::left_join(., precipdf)


#......................................................................................................
# Precipitation and Temperature Considered as Means
#......................................................................................................
# months of study period
studyperiod <- levels(factor(wthrnd.mnth$hvyrmnth_dtmnth))

precipdf <- precipdf %>% 
  dplyr::filter(hvyrmnth_dtmnth %in% studyperiod)

precipstack <- raster::stack(precipdf$precipraster)
precipstack.mean <- raster::calc(precipstack, mean, na.rm = T)

tempdf <- tempdf %>% 
  dplyr::filter(hvyrmnth_dtmnth %in% studyperiod)

tempstack <- raster::stack(tempdf$tempraster)
tempstack.mean <- raster::calc(tempstack, mean, na.rm = T)

wthrnd.mean <- dt[,c("hv001", "geometry", "urban_rura")] %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd.mean <- wthrnd.mean[!duplicated(wthrnd.mean$hv001),]


# Drop in a for loop again to account for DHS buffering
wthrnd.mean$precip_mean_cont_clst <- NA
wthrnd.mean$temp_mean_cont_clst <- NA

for(i in 1:nrow(wthrnd.mean)){
  # precip
  wthrnd.mean$precip_mean_cont_clst[i] <- 
    raster::extract(x = precipstack.mean, # this doesn't change this time 
                    y = sf::as_Spatial(wthrnd.mean$geometry[i]),
                    buffer = wthrnd.mean$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
  
  # temp
  wthrnd.mean$temp_mean_cont_clst[i] <- 
    raster::extract(x = tempstack.mean, # this doesn't change this time 
                    y = sf::as_Spatial(wthrnd.mean$geometry[i]),
                    buffer = wthrnd.mean$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
  
}

wthrnd.mean <- wthrnd.mean %>% 
  dplyr::select(c("hv001", "precip_mean_cont_clst", "temp_mean_cont_clst"))
sf::st_geometry(wthrnd.mean) <- NULL

#............................................................................
# OUT
#............................................................................
saveRDS(object = precipstack.mean, 
        file = "data/derived_data/vividepi_precip_study_period_effsurface.rds")


saveRDS(object = wthrnd.mean, 
        file = "data/derived_data/vividep_weather_recoded_mean.rds")






