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
sp::proj4string(caf) <- "+init=epsg:4326"
# need mask 
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")
sf::st_crs(DRC)

#......................................................................................................
# First make Temperature raster for air temperature from Garske's Estimating Air Temperature and Its Influence on Malaria Transmission across Africa
# PMCID: PMC3577915
#......................................................................................................
#......................
# read in previously published data and follow publicly available scripting 
#......................
# parameters used in the Fourier Transform:
nyears = 4 # used 4 years of data to calculate the Fourier Transform
ppyear = 64 # used a time step of 64 points per year in the input data to calculate the Fourier Transforms.

# reading in the Fourier Transform data file:
FT = read.table("data/raw_data/prev_pub/paper_airTemperature/FT_2003-2006_64ppyear_ERAday_0_1deg_255721.txt",header=T)

# specifying the time points throughout the year
# for which to calculate the time series.
# here: daily time step
timepoints = (1:365)/365

# calculating the constant term (=annual mean) for each location and time point:
TS0 = matrix(FT$H0/(nyears*ppyear),nrow=nrow(FT),ncol=length(timepoints))
# calculating the time series based on
# constant term H0 and annual mode H1:
TS1 = TS0 + 1/(nyears*ppyear)*(outer(FT$ReH1,cos(2*pi*timepoints*1))+outer(FT$ImH1,sin(2*pi*timepoints*1)))
# calculating the time series based on
# constant term H0, annual mode H1 and biannual mode H2
TS2 = TS1 + 1/(nyears*ppyear)*(outer(FT$ReH2,cos(2*pi*timepoints*2))+outer(FT$ImH2,sin(2*pi*timepoints*2)))

#......................
# will use annual mode H1 and biannual mode H2
# create africa rasters
#......................
aftemp <- cbind.data.frame(FT[,c("longitude", "latitude")], TS2)
dayrstr_list <- list()
for(i in 3:ncol(aftemp)) {
  dayrstr <- raster::rasterFromXYZ(xyz = aftemp[,c(1, 2, i)],
                                   crs = "+init=epsg:4326")
  dayrstr_list <- append(dayrstr_list, dayrstr)
}
# stack
daytemprstr_stack <- raster::stack(dayrstr_list)
names(daytemprstr_stack) <- paste0("X", 1:365)
# crop and mask
daytemprstr_stack <- raster::crop(daytemprstr_stack, caf)
daytemprstr_stack <- raster::mask(daytemprstr_stack, DRC)


#......................................................................................................
# Get Precipitation Rasters from (CHIRPS)
#......................................................................................................
readRasterBB.precip <- function(rstfile, sp = sp, caf = caf){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = caf)
  ret <- raster::mask(x = ret, mask = sp)
  
  vals <- raster::values(ret) 
  vals <- ifelse(vals == -9999, NA, vals) # improper values
  raster::values(ret) <- vals
  
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+init=epsg:4326")) 
  
  return(ret)
  
}
#......................
# files and read in
#......................
precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T, 
                     pattern = ".tif")
precipfrst <- lapply(precip, readRasterBB.precip, sp = DRC, caf = caf)

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

#......................................................................................................
# Precipitation and Temperature Considered as Means
#......................................................................................................
#......................
# Get study months to subset rasters
#   and tidy up for means later
#......................
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# read in GE as import
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")
# drop observations with missing or excluded geospatial data 
dt <- dt  %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>%  # drop observations with missing geospatial data 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>%
  dplyr::filter(hv102 == 1) %>% # subset to de-jure https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm
  dplyr::filter(hiv05 != 0) # sampling weights 0

# sanity check
sf::st_crs(dt)
identicalCRS(sf::as_Spatial(dt), sf::as_Spatial(ge))
identicalCRS(sf::as_Spatial(dt), caf)
# liftover to conform with rgdal updates http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
dt <- sp::spTransform(sf::as_Spatial(dt), CRSobj = sp::CRS("+init=epsg:4326"))
identicalCRS(dt, caf)
# back to tidy 
dt <- sf::st_as_sf(dt)

#......................
# subset to correct "dates"
#......................
studyperiod <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
  dplyr::select(c("hvdate_dtdmy", "hvyrmnth_dtmnth")) %>% 
  dplyr::filter(!duplicated(.))
sf::st_geometry(studyperiod) <- NULL

# months of study period
mtnhs <- unique(studyperiod$hvyrmnth_dtmnth)

# subset precipdf 
precipdf <- precipdf %>% 
  dplyr::filter(hvyrmnth_dtmnth %in% mtnhs)

# make precip df stack mean
precipstack <- raster::stack(precipdf$precipraster)
precipstack.mean <- raster::calc(precipstack, mean, na.rm = T)

# for temp
studyperiod <- studyperiod %>% 
  dplyr::mutate(day = ifelse(lubridate::year(hvdate_dtdmy) == 2014,
                             hvdate_dtdmy - lubridate::ymd("20140101") + 1,
                             hvdate_dtdmy - lubridate::ymd("20130101") + 1))
daytemprstr_stack <- daytemprstr_stack[[paste0("X", unique(sort(studyperiod$day)))]]

# make temp df stack mean
daytemprstr_stack.mean <- raster::calc(daytemprstr_stack, mean, na.rm = T)

#......................
# now lets get cluster level means
#......................
wthrnd.mean <- ge[,c("hv001", "geometry", "urban_rura")] %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10000, 2000))


# Drop in a for loop again to account for DHS buffering
# note the 0.05 degree resolution for precip is approximately 5.5km, so the buffer
# note, the 0.1 degree resolution for temp is approx 11km, so exceeds buffer but will do for consistency
# for urbanicity shouldn't be doing anything... but to be consistent with 
# DHS "The Geospatial Covariate Datasets Manual", we will do it
sf::st_crs(precipstack.mean)
sf::st_crs(daytemprstr_stack.mean)
sf::st_crs(wthrnd.mean)
identicalCRS(daytemprstr_stack.mean, caf)
identicalCRS(daytemprstr_stack.mean, precipstack.mean)
identicalCRS(daytemprstr_stack.mean, sf::as_Spatial(wthrnd.mean))

# now run for loop
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
    raster::extract(x = daytemprstr_stack.mean, # this doesn't change this time 
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
    raster::extract(x = daytemprstr_stack.mean, # this doesn't change this time 
                    y = sf::as_Spatial(wthrnd.mean$geometry[i]),
                    buffer = 12000,
                    fun = mean,
                    na.rm = T, 
                    sp = F
    )
}


#......................
# check outputs
#......................
wthrnd.mean %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2]) %>% 
  ggplot() + 
  geom_sf(data = sf::st_as_sf(DRC)) +
  ggspatial::layer_spatial(data = precipstack.mean,
                           aes(fill = stat(band1)),
                           alpha = 0.9,
                           na.rm = T) +
  geom_point(aes(x = longnum, y = latnum, color = precip_mean_cont_clst)) +
  scale_fill_viridis_c(option="plasma", direction = 1) +
  scale_color_viridis_c(option="viridis") + 
  ggtitle("Precip")


wthrnd.mean %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2]) %>% 
  ggplot() + 
  geom_sf(data = sf::st_as_sf(DRC)) +
  ggspatial::layer_spatial(data = daytemprstr_stack.mean,
                           aes(fill = stat(band1)),
                           alpha = 0.9,
                           na.rm = T) +
  geom_point(aes(x = longnum, y = latnum, color = temp_mean_cont_clst)) +
  scale_fill_viridis_c(option="plasma", direction = 1) +
  scale_color_viridis_c(option="viridis") + 
  ggtitle("Temp")


#............................................................................
# OUT
#............................................................................
wthrnd.mean <- wthrnd.mean %>% 
  dplyr::select(c("hv001", "precip_mean_cont_clst", "temp_mean_cont_clst"))
sf::st_geometry(wthrnd.mean) <- NULL


saveRDS(object = precipstack.mean, 
        file = "data/derived_data/vividepi_precip_study_period_effsurface.rds")

saveRDS(object = daytemprstr_stack.mean, 
        file = "data/derived_data/vividepi_temperature_study_period_effsurface.rds")

saveRDS(object = wthrnd.mean, 
        file = "data/derived_data/vividep_weather_recoded_mean.rds")









