





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







