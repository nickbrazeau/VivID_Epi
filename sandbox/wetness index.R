
#---------------------------------------------------------------------------------
# Wetness Index
#---------------------------------------------------------------------------------
env <- rsaga.env()
env$workspace <- "~/Documents/GitHub/VivID_Epi/data/raw_data/saga_analysis/"

drcdem <- elevatr::get_elev_raster(sf::as_Spatial(bb), z=6)
raster::writeRaster(drcdem, filename = "data/raw_data/saga_analysis/drcdem.sgrd", 
                    format = "SAGA", overwrite = T)

rsaga.wetness.index(in.dem = "drcdem", 
                    out.wetness.index = "drcdem_wetnessindex", 
                    env = env)
sagawetness <- raster(x = "data/raw_data/saga_analysis/drcdem.sdat") # read in wetness

# smooth and liftover raster for ggplot 
sagawetness.pts  <-  raster::rasterToPoints(sagawetness)
sagawetness.pts.df <-  data.frame(lon = sagawetness.pts[,1], 
                                  lat = sagawetness.pts[,2], 
                                  wetness = sagawetness.pts[,3])
ggplot() + 
  geom_raster(data = sagawetness.pts.df, 
              aes(lon, lat, fill = wetness), alpha = 0.5) +
  geom_sf(data=DRCprov) +
  scale_fill_gradient2("Wetness", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") 

