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
# Set up Lagged Month for Study Collection time
#......................................................................................................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")


# drop clusters with missing geospatial data 
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 

#.............
# dates
#.............
dt <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")))


# NOTE, some clusters have survey start and end dates that are in two months 
# (eg boundaries aren't clean/coinciding with a month. Naturally). Given
# grouping by month, need to assign a clusters "month" on the majority of days 
# that were spent surveying that clusters

# clusters without clean boundaries
clst_mnth_bounds <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
  dplyr::mutate(mnth = lubridate::month(hvdate_dtdmy)) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(moremnths = length(unique(mnth))) %>% 
  dplyr::filter(moremnths > 1)

clst_mnth_bounds.assign <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
  dplyr::filter(hv001 %in% clst_mnth_bounds$hv001) %>% 
  dplyr::mutate(hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
  dplyr::group_by(hv001, hvyrmnth_dtmnth) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n == max(n)) %>% 
  dplyr::select(-c("n"))

sf::st_geometry(clst_mnth_bounds.assign) <- NULL

# first join assigned months
# and then make date for unambigious months
dt <- dt %>% 
  dplyr::left_join(x=., y = clst_mnth_bounds.assign, by = "hv001") %>% 
  dplyr::mutate(hvyrmnth_dtmnth = ifelse(is.na(hvyrmnth_dtmnth),
                                         paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"),
                                         hvyrmnth_dtmnth))

dates <- readr::read_csv("internal_datamap_files/pr_date_liftover.csv")
dt <- dt %>% 
  dplyr::left_join(x=., y=dates, by = "hvyrmnth_dtmnth") %>% 
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))

xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)


# final dataframe for lagged months
clst_mnths.lag <- dt %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag"))
sf::st_geometry(clst_mnths.lag) <- NULL  
clst_mnths.lag <- clst_mnths.lag %>% 
  dplyr::filter(!duplicated(.))

# months of study period
studyperiod <- levels(factor(dt$hvyrmnth_dtmnth))

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
                hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                precipraster = precipfrst) %>% 
  dplyr::select(c("hvyrmnth_dtmnth_lag", "precipraster"))



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
                hvyrmnth_dtmnth_lag = factor(paste(year, month, sep = "-")),
                tempraster = tempfrst
  ) %>% 
  dplyr::select(c("hvyrmnth_dtmnth_lag", "tempraster"))




#......................................................................................................
# Precipitation and Temperature Considered by Month 
#......................................................................................................

wthrnd.mnth <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd.mnth <- wthrnd.mnth[!duplicated(wthrnd.mnth$hv001),]

wthrnd.mnth <- wthrnd.mnth %>% 
  dplyr::left_join(., tempdf) %>% 
  dplyr::left_join(., precipdf)


# Drop in a for loop to acount for dhs buffering
wthrnd.mnth$precip_lag_cont_clst <- NA
wthrnd.mnth$temp_lag_cont_clst <- NA

for(i in 1:nrow(wthrnd.mnth)){
  # precip
  wthrnd.mnth$precip_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd.mnth$precipraster[[i]],
                    y = sf::as_Spatial(wthrnd.mnth$geometry[i]),
                    buffer = wthrnd.mnth$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
  
  # temp
  wthrnd.mnth$temp_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd.mnth$tempraster[[i]],
                    y = sf::as_Spatial(wthrnd.mnth$geometry[i]),
                    buffer = wthrnd.mnth$buffer[i],
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
  
}

wthrnd.mnth <- wthrnd.mnth %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst", "temp_lag_cont_clst")) %>% 
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
sf::st_geometry(wthrnd.mnth) <- NULL


wthrnd.mnth <- dplyr::inner_join(clst_mnths.lag, wthrnd.mnth, by = c("hv001", "hvyrmnth_dtmnth_lag")) # dominant month needs to win




#......................................................................................................
# Precipitation and Temperature Considered as Means
#......................................................................................................
precipdf <- precipdf %>% 
  dplyr::filter(hvyrmnth_dtmnth_lag %in% studyperiod)

precipstack <- raster::stack(precipdf$precipraster)
precipstack.mean <- raster::calc(precipstack, mean, na.rm = T)

tempdf <- tempdf %>% 
  dplyr::filter(hvyrmnth_dtmnth_lag %in% studyperiod)

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
saveRDS(object = wthrnd.mnth, 
        file = "data/derived_data/vividep_weather_recoded_monthly.rds")

saveRDS(object = wthrnd.mean, 
        file = "data/derived_data/vividep_weather_recoded_mean.rds")






