library(tidyverse)
library(sf)

#.............................................................
# Precipitation (CHRIPS) 
#.............................................................

# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+proj=longlat +datum=WGS84 +no_defs"

# precip data
precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T, 
                     pattern = ".tif")
precipfrst <- lapply(precip, readRasterBB, bb = caf)

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


#.............................................................
# Temperature (MODIS/LAADS) 
#.............................................................
# FILENAME DECONVOLUTION: An example of an existing file is below. 
# https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MOD02QKM/2007/018/MOD02QKM.A2007018.0105.006.2014227230926.hdf
# This path should return a MODIS Terra quarter kilometer (250 m) top of atmosphere reflectance product for year 2007, day-of-year 018 (i.e. January 18), from collection 6.
# Note to change values back to celsius: https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius

# for temperature, need to mask water sources which have too low temp readings
load("data/derived_data/hotosm_waterways.RDA")
oceans <- sf::st_read("data/map_bases/ne_10m_ocean/ne_10m_ocean.shp") %>% 
  sf::st_crop(caf)

tempfiles <- list.files(path = "data/raw_data/weather_data/LAADS_NASA/", full.names = T,
                        pattern = "_Night_CMG.tif")
tempfrst <- lapply(tempfiles, readRasterBB, bb = caf)

# rescale values to celsius
tempfrst <- lapply(tempfrst, function(x){
  # mask water
  x <- raster::mask(x, sf::as_Spatial(wtrlns))
  x <- raster::mask(x, sf::as_Spatial(wtrply))
  x <- raster::mask(x, sf::as_Spatial(oceans))
  
  
  vals <- raster::values(x) 
  vals <- ifelse(vals <= 7500, NA, vals) # improper values
  vals <- (vals * 0.02) - 273.15
  raster::values(x) <- vals
  return(x)
})

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


wthrnd <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd <- wthrnd[!duplicated(wthrnd$hv001),]

wthrnd <- wthrnd %>% 
  dplyr::left_join(., tempdf) %>% 
  dplyr::left_join(., precipdf)

# Drop in a for loop
wthrnd$precip_lag_cont_clst <- NA
wthrnd$temp_lag_cont_clst <- NA

for(i in 1:nrow(wthrnd)){
  # precip
  wthrnd$precip_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd$precipraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )
  
  # temp
  wthrnd$temp_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd$tempraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )
  
}

wthrnd <- wthrnd %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst", "temp_lag_cont_clst")) %>% 
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
sf::st_geometry(wthrnd) <- NULL

#........................
# OUT
#........................
saveRDS(object = wthrnd, 
        file = "data/derived_data/vividep_weather_recoded.rds")






