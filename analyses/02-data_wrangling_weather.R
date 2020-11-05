#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle weather data 
# that is around the time of our study period for the CD2013
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(raster)
library(sp)
library(sf)

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"



#..............................
# Housekeeping
#..............................
readRasterBB.precip <- function(rstfile, sp = sp, caf = caf){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = caf)
  ret <- raster::mask(x = ret, mask = sp)
  
  vals <- raster::values(ret) 
  vals <- ifelse(vals == -9999, NA, vals) # improper values
  raster::values(ret) <- vals
  
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs +units=m")) # want units to be m
  
  return(ret)
  
}

readRasterBB.temp <- function(rstfile, sp = sp, caf = caf){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = caf)
  ret <- raster::mask(x = ret, mask = sp)
  
  vals <- raster::values(ret) 
  vals <- ifelse(vals <= 7500, NA, vals) # improper values
  vals <- (vals * 0.02) - 273.15
  raster::values(ret) <- vals
  
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs +units=m")) # want units to be m
  return(ret)
}


# create mask 
DRCprov <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")

#......................................................................................................
# Precipitation (CHRIPS) and Temperature (MODIS/LAADS) Read In Data
#......................................................................................................

precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T, 
                     pattern = ".tif")
precipfrst <- lapply(precip, readRasterBB.precip, sp = DRCprov, caf = caf)

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
tempfrst <- lapply(tempfiles, readRasterBB.temp, sp = DRCprov, caf = caf)

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
# Precipitation and Temperature Considered as Means
#......................................................................................................
# Get study months
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# drop observations with missing geospatial data 
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 

wthrnd.mnth <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth"))


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
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))
wthrnd.mean <- wthrnd.mean[!duplicated(wthrnd.mean$hv001),]


# Drop in a for loop again to account for DHS buffering
# note the 0.05 degree resolution is approximately 6km, so the buffer
# for urbanicity shouldn't be doing anything... but to be consistent with 
# DHS "The Geospatial Covariate Datasets Manual", we will do it
sf::st_crs(precipstack.mean)
sf::st_crs(tempstack.mean)
sf::st_crs(dt)


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


#..............................................................
# Interpolate missing boundaries
#..............................................................
# sanity check for missing
missingclst <- wthrnd.mean %>% 
  dplyr::filter(c( is.na(precip_mean_cont_clst) | is.na(temp_mean_cont_clst) )) %>% 
  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])

ggplot() +
  geom_sf(data = wthrnd.mean) +
  geom_point(data = missingclst, aes(x=long, y=lat), color = "red") +
  ggrepel::geom_label_repel(data = missingclst, aes(x=long, y=lat, label = hv001)) 

# all boundaries and small buffer (2km is barely outside raster cell often), will increase to 6km
precipmissing <- which(is.na(wthrnd.mean$precip_mean_cont_clst))
for (i in precipmissing) {
  # precip
  wthrnd.mean$precip_mean_cont_clst[i] <- 
    raster::extract(x = precipstack.mean, # this doesn't change this time 
                    y = sf::as_Spatial(wthrnd.mean$geometry[i]),
                    buffer = 6000,
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}

tempmissing <- which(is.na(wthrnd.mean$temp_mean_cont_clst))
for (i in tempmissing) {
  # temp
  wthrnd.mean$temp_mean_cont_clst[i] <- 
    raster::extract(x = tempstack.mean, # this doesn't change this time 
                    y = sf::as_Spatial(wthrnd.mean$geometry[i]),
                    buffer = 6000,
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

saveRDS(object = tempstack.mean, 
        file = "data/derived_data/vividepi_temperature_study_period_effsurface.rds")

saveRDS(object = wthrnd.mean, 
        file = "data/derived_data/vividep_weather_recoded_mean.rds")






