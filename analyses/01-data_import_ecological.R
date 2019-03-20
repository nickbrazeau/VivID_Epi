#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import ECOLOGICAL Variables from the web that has to do with geospatial and climate data
# Will then merge to dhs scrape in the 01-data_import_epi file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
devtools::install_github("ropensci/osmdata")
library(osmdata)
devtools::install_github("malaria-atlas-project/malariaAtlas")
library(malariaAtlas)
devtools::install_github("OJWatson/magenta")
library(magenta)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")

#---------------------------------------------------------------------------------
# set up DRC borders
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata

bb <- getbb("Democratic Republic of the Congo", featuretype = "country")
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')

#---------------------------------------------------------------------------------
# Pull Down Great Ape Territories from IUC (minus Pongo)
#---------------------------------------------------------------------------------
primate <- sf::read_sf("data/redlist_species_data_primate/data_0.shp")
ape <- primate[grepl("pan paniscus|pan troglodytes|gorilla", tolower(primate$BINOMIAL)), ] # pan trog, pan panisus, gorilla sp


#---------------------------------------------------------------------------------
# Seasonality from Imperial/Carins
#---------------------------------------------------------------------------------

# From OJ (via slack on 1/10/2019)
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056487
# https://www.nature.com/articles/ncomms1879#supplementary-information
# https://betterexplained.com/articles/an-interactive-guide-to-the-fourier-transform/
admin_units_seasonal <- readRDS("data/imperial_share/admin1_seasonality_for_nick_from_OJ.rds")
drc_admin_units_seasonal <- admin_units_seasonal %>% 
  filter(COUNTRY_ID == "COD")

### NEED TO FINISH -- ERROR WITH THETA_C still
source("~/Documents/GitHub/VivID_Epi/R/00-functions_seasonality.R")




#---------------------------------------------------------------------------------
# Rainfall Data from OJ for admin1 in DHS
#---------------------------------------------------------------------------------
# Daily DRC's real rainfall data from CHIRPs for 2013-2017

rain <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/imperial_share/drc_rainfall_2013-2017.rds")

 
#---------------------------------------------------------------------------------
# MAP project for R 
#---------------------------------------------------------------------------------
# from their example code chunk
# download raster of travel time to cities (Weiss et al 2018) for study area & visualise this
TravelToCities <- malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                                  extent = bb)
TravelToCities <- log(TravelToCities+0.1) # MAP project decided 0.1 for offest. going to keep it



#---------------------------------------------------------------------------------
# Link to GE
#---------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------
# QUICK Interpolation for Temperature and Precip
#----------------------------------------------------------------------------------------------------
# https://mgimond.github.io/Spatial/interpolation-in-r.html
tempprecip <- dt %>% 
  group_by(hv001) %>% 
  summarise(clsttemp = mean(mean_temperature_2015_cont),
            clstsdtemp = sd(mean_temperature_2015_cont),
            clstrain = mean(rainfall_2015_cont),
            clstsdrain = sd(rainfall_2015_cont)) %>% 
  left_join(., y=clustgeom, by = "hv001")

sum(apply(tempprecip[,c(3,5)], 1, function(x){return(x > 0)}), na.rm = T)


# imputing for mean
tempprecip$clsttemp[is.na(tempprecip$clsttemp)] <- mean(tempprecip$clsttemp, na.rm=T)

#---------------------------------------
# Temperature
#---------------------------------------

#.....................
# make raster and boundaries
#....................
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long", "lat")
pos <- sf::as_Spatial(sf::st_as_sf(polybb))

st_sf(bb = 1:2, geom = st_sfc( st_point(bb[,1]), st_point(bb[,2]) ) )


grd              <- as.data.frame(spsample(grid.pred, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(pos)

temp.idw <- gstat::idw(tempprecip$clsttemp ~ 1, pos, 
                       newdata = grd, idp = 2.0)
# Convert to raster object then clip to DRC
temp.r       <- raster(temp.idw)
temp.r.m     <- mask(temp.r, as_Spatial(DRCprov))

# First smooth and then make the temperature raster plot
temp.r.m <- focal(temp.r.m, w=matrix(1,
                                     nrow=5,
                                     ncol=5), mean)
temp.r.m.pts  <-  raster::rasterToPoints(temp.r.m)
temp.r.m.pts.df <-  data.frame(lon = temp.r.m.pts[,1], 
                               lat = temp.r.m.pts[,2], 
                               temp = temp.r.m.pts[,3])
temp.r.m.pts.m.plot <- ggplot() + 
  geom_raster(data = temp.r.m.pts.df, aes(lon, lat, fill = temp), alpha = 0.8) +
  scale_fill_gradient2("Temperature", low = "#313695", mid = "#ffffbf", high = "#a50026",
                       midpoint = 21) +
  prettybasemap_nodrc 



jpeg("~/Documents/GitHub/VivID_Epi/figures/07-temperature-idw.jpg", 
     height = 8, width=8, units = "in", res=400)
plot(temp.r.m.pts.m.plot)
graphics.off()


#---------------------------------------
# Precipation
#---------------------------------------
precip.idw <- gstat::idw(tempprecip$clstrain ~ 1, pos, 
                         newdata = grd, idp = 2.0)
# Convert to raster object then clip to DRC
precip.r       <- raster(precip.idw)
precip.r.m     <- mask(precip.r, as_Spatial(DRCprov))

# First smooth and then make the temperature raster plot
precip.r.m <- focal(precip.r.m, w=matrix(1,
                                         nrow=5,
                                         ncol=5), mean)
precip.r.m.pts  <-  raster::rasterToPoints(precip.r.m)
precip.r.m.pts.df <-  data.frame(lon = precip.r.m.pts[,1], 
                                 lat = precip.r.m.pts[,2], 
                                 precip = precip.r.m.pts[,3])
precip.r.m.pts.m.plot <- ggplot() + 
  geom_raster(data = precip.r.m.pts.df, aes(lon, lat, fill = precip), alpha = 0.8) +
  scale_fill_gradient2("Precipitation", low = "#313695", mid = "#ffffbf", high = "#a50026",
                       midpoint = 1520)+
  prettybasemap_nodrc 



jpeg("~/Documents/GitHub/VivID_Epi/figures/07-precip-idw.jpg", 
     height = 8, width=8, units = "in", res=400)
plot(precip.r.m.pts.m.plot)
graphics.off()







