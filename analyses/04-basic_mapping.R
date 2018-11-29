# add in GPS points
cd2013gps <- tibble(cluster = cd2013dhsge$DHSCLUST, lat = cd2013dhsge$LATNUM, long = cd2013dhsge$LONGNUM)
# Clusters with Lat and Long of 0,0 were not able to be identified and should have coordinates set to NA
cd2013gps <- cd2013gps %>%
  dplyr::mutate(long = ifelse(long == 0 & lat == 0, NA, long)) %>%
  dplyr::mutate(lat = ifelse(long == 0 & lat == 0, NA, lat))