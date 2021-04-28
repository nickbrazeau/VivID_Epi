#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle 
# the covars needed to fit the prov carbayes models
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
library(raster)
library(sf)
source("R/00-functions_basic.R")
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")

#......................
# get and wrangle hlthdist
#......................
# items for mask
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")

# read in
hlthdist <- raster::raster("data/raw_data/hlthdist/2020_walking_only_travel_time_to_healthcare.geotiff")
# sanity check
sf::st_crs(hlthdist)
identicalCRS(hlthdist, caf)
identicalCRS(hlthdist, sf::as_Spatial(DRCprov))

# crop for speed
hlthdist <- raster::crop(x = hlthdist, y = caf)
# create mask 
hlthdist.drc <- raster::mask(x = hlthdist, mask = DRCprov)


#-------------------------------------------------------------------------
# Aggregate Covariates
# means within each adm1
#-------------------------------------------------------------------------
precip <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds")
hlthdist.drc 

extract_agg_raster_polygon <- function(rstrlyr, plygn){
  vals <- raster::extract(x = rstrlyr, y = sf::as_Spatial(plygn),
                          fun = mean,
                          na.rm = T, 
                          sp = F
  )
  return(as.vector(vals))
  
}

# split admins
adm1 <- DRCprov %>% 
  dplyr::select(c("adm1name", "geometry")) %>% 
  dplyr::filter(!duplicated(.))

# sanity
sp::identicalCRS(sf::as_Spatial(adm1), caf)
sp::identicalCRS(sf::as_Spatial(adm1), precip)
sp::identicalCRS(sf::as_Spatial(adm1), hlthdist.drc)

# store names
pvcovar <- adm1 %>% 
  dplyr::select(c("adm1name"))

# split into list
adm1 <- split(adm1, 1:nrow(adm1))


# get covars 
pvcovar$precip <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = precip) )
pvcovar$hlthdist <- unlist( lapply(adm1, extract_agg_raster_polygon, rstrlyr = hlthdist.drc) )

#...................... 
# sanity checks 
#......................
# quick viz
hist(pvcovar$precip)
hist(pvcovar$hlthdist)

#  precip 
pa <- ggplot() +
  geom_sf(data = pvcovar, aes(fill = precip)) + 
  scale_fill_viridis_c("Precip")
pb <- ggplot() + 
  ggspatial::layer_spatial(data = precip, aes(fill = stat(band1)),
                           alpha = 0.8) +
  geom_sf(data = DRCprov, color = "#000000", fill = NA) +
  scale_fill_viridis_c("Precip", na.value = NA) 

cowplot::plot_grid(pa, pb)



# health distance 
pa <- ggplot() +
  geom_sf(data = pvcovar, aes(fill = hlthdist)) + 
  scale_fill_viridis_c("hlthdist")
pb <- ggplot() + 
  ggspatial::layer_spatial(data = hlthdist.drc, aes(fill = stat(band1)),
                           alpha = 0.8) +
  geom_sf(data = DRCprov, color = "#000000", fill = NA) +
  scale_fill_viridis_c("Health Distance", na.value = NA) 

jpeg("~/Desktop/temp.jpg", width = 11, height = 8, res = 250, units = "in")
cowplot::plot_grid(pa, pb)
graphics.off()



#...................... 
# scale for better model fitting
#......................
pvcovar$precip_scale <- my.scale(pvcovar$precip)
pvcovar$hlthdist_scale <- my.scale(pvcovar$hlthdist)
sf::st_geometry(pvcovar) <- NULL




#..............................................................
# save out
#..............................................................
saveRDS(pvcovar, file = "data/derived_data/vividepi_prov_covars_bayesian_fit.RDS")
