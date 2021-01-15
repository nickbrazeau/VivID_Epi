#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get spatial predictions from 
# spatial models produced from prevmap
#----------------------------------------------------------------------------------------------------
source("R/00-functions_maps.R") 
source("R/00-functions_basic.R")
library(tidyverse)
library(drake)
library(raster)
library(PrevMap)
set.seed(48, "L'Ecuyer")

#......................
# pieces need to manipulate spatial preds
#......................
# create bounding box of Central Africa for Speed
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"

# create mask 
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds")


#...............................
# sample coordinates for prediction surface 
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.coords <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred.coords) <- c("longnum","latnum")
grid.pred.coords.df <- as.data.frame(grid.pred.coords)
# convert to sf tidy  
grid.pred.coords <- sf::st_as_sf(grid.pred.coords.df, coords = c("longnum", "latnum"),
                                 crs = "+init=epsg:4326")

#...............................
# extract covariate information for prediction surface 
#...............................
# Raster surfaces for risk factors
riskvars <- c("precip_mean_cont_scale_clst", 
              "cropprop_cont_scale_clst", "hlthdist_cont_scale_clst")

# take mean aggregate up and then scale like covar
precipraster <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 
precipraster <- raster::aggregate(precipraster, fact = 2, fun = mean)
# transform like fitting
precipraster <- raster::scale(precipraster, center = TRUE, scale = TRUE)

# take mean here since we converted to binary
cropraster <- readRDS("data/derived_data/vividepi_cropland_surface.rds")
cropraster <- raster::aggregate(cropraster, fact = 36, fun = mean)
# transform like fitting
values(cropraster) <- my.scale(logit(values(cropraster), tol = 1e-3))


# here we need to take healthcare walk distance
# aggregate up and then do log transform so it more closelymatches how we manipulate before
htlhdistraster <- raster::raster("data/raw_data/hlthdist/2020_walking_only_travel_time_to_healthcare.geotiff")
sp::identicalCRS(htlhdistraster, caf)
sp::identicalCRS(htlhdistraster, sf::as_Spatial(DRC))
htlhdistraster <- raster::crop(htlhdistraster, caf)
htlhdistraster <- raster::mask(htlhdistraster, DRC)
htlhdistraster <- raster::aggregate(htlhdistraster, fact = 12, fun = mean)
values(htlhdistraster) <- my.scale(log(values(htlhdistraster) + 0.5))

# stack 
predcovars <- raster::stack(precipraster, cropraster, htlhdistraster)
names(predcovars) <- riskvars

pred.df <- raster::extract(
  x = predcovars,
  y = sf::as_Spatial(grid.pred.coords),
  buffer = 10000, # max dhs offset
  fun = mean,
  na.rm = T,
  sp = F)

pred.df <- as.data.frame(pred.df)
pred.df <- cbind.data.frame(grid.pred.coords.df, pred.df)

#......................
# sanity checks
#......................
p1 <- ggplot() +
  ggspatial::layer_spatial(data = precipraster, aes(fill = stat(band1))) +
  scale_fill_viridis_c("", na.value = NA) +
  theme_void()
p2 <- ggplot() +
  geom_point(data = pred.df, aes(x = longnum, y = latnum,
                                 color = precip_mean_cont_scale_clst))  +
  scale_color_distiller("", type = "div", palette = "RdYlBu", na.value = NA) +
  theme_void()
cowplot::plot_grid(p1, p2, nrow = 1)

p1 <- ggplot() +
  ggspatial::layer_spatial(data = cropraster, aes(fill = stat(band1))) +
  scale_fill_viridis_c("", na.value = NA) +
  theme_void()
p2 <- ggplot() +
  geom_point(data = pred.df, aes(x = longnum, y = latnum,
                                 color = cropprop_cont_scale_clst))  +
  scale_color_distiller("", type = "div", palette = "RdYlBu", na.value = NA) +
  theme_void()
cowplot::plot_grid(p1, p2, nrow = 1)

p1 <- ggplot() +
  ggspatial::layer_spatial(data = htlhdistraster, aes(fill = stat(band1))) +
  scale_fill_viridis_c("", na.value = NA) +
  theme_void()
p2 <- ggplot() +
  geom_point(data = pred.df, aes(x = longnum, y = latnum,
                                 color = hlthdist_cont_scale_clst))  +
  scale_color_distiller("", type = "div", palette = "RdYlBu", na.value = NA) +
  theme_void()
cowplot::plot_grid(p1, p2, nrow = 1)

# look at hlth dist which from MAP was very pointed
p1 <- ggplot() +
  ggspatial::layer_spatial(data = htlhdistraster_raw, aes(fill = stat(band1))) +
  scale_fill_viridis_c("", na.value = NA) +
  theme_void()
p2 <- ggplot() +
  ggspatial::layer_spatial(data = htlhdistraster, aes(fill = stat(band1))) +
  scale_fill_viridis_c("", na.value = NA) +
  theme_void()
cowplot::plot_grid(p1, p2, nrow = 1)


#...............................
# Bounds so we only do interpolation
# (i.e. coerce extrapolation values)
#...............................
# precip
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
precip <- dt %>% 
  dplyr::select("hv001", "precip_mean_cont_scale_clst") %>% 
  dplyr::filter(!duplicated(.))
max.precip <- max(precip$precip_mean_cont_scale_clst)
min.precip <- min(precip$precip_mean_cont_scale_clst)

pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst > max.precip] <- max.precip
pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst < min.precip] <- min.precip


# bring in cropland
crop <- readRDS("data/derived_data/vividepi_cropland_propmeans.rds") 
max.crop <- max(crop$cropprop_cont_scale_clst)
min.crop <- min(crop$cropprop_cont_scale_clst)
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst > max.crop] <- max.crop
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst < min.crop] <- min.crop


# bring in health care distance/access clst means
acc <- readRDS("data/derived_data/vividepi_hlthdist_clstmeans.rds") 
max.acc <- max(acc$hlthdist_cont_scale_clst)
min.acc <- min(acc$hlthdist_cont_scale_clst)
pred.df$hlthdist_cont_scale_clst[pred.df$hlthdist_cont_scale_clst > max.acc] <- max.acc
pred.df$hlthdist_cont_scale_clst[pred.df$hlthdist_cont_scale_clst < min.acc] <- min.acc


#..............................................................
# Convert to raster for sampling
#..............................................................
covar.rstr.pred <- raster::rasterFromXYZ(pred.df, 
                                         res = c(0.1, 0.1),
                                         crs = "+init=epsg:4326")

#......................
# save out
#......................
# save this out for spatial preds covar
saveRDS(covar.rstr.pred, file = "data/derived_data/vividepi_spatial_covar_feature_engineer.rds")
