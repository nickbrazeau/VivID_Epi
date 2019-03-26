

# From the MAP Project
# download raster of travel time to cities (Weiss et al 2018) for study area & visualise this
bb <- osmdata::getbb("Democratic Republic of the Congo", 
                     featuretype = "country")

TravelToCities <- malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                                          extent = bb)
TravelToCities <- raster::projectRaster(from = TravelToCities, to = TravelToCities,
                                        crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m


clstsrch <- dt %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2)) %>%  # after DHS GC manual
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(buffer = mean(buffer))


clstsrch$travel_mean <- raster::extract(x = TravelToCities, 
                                        y = sf::as_Spatial( clstsrch$geometry ), 
                                        buffer = clstsrch$buffer,
                                        fun = mean,
                                        sp = F) 

clstsrch <- clstsrch %>% 
  dplyr::select(-c("geometry")) %>% 
  as.data.frame(.) # for easier merge below