#----------------------------------------------------------------------------------------------------
# Purpose of this script is partition the DRC
# by geographic coordinates for spatial cross-validation
#----------------------------------------------------------------------------------------------------

library(tidyverse)
library(sf)


DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds"))

coords <- ge %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::select(c("dhsclust", "longnum", "latnum")) %>% 
  dplyr::filter(!duplicated(.)) %>% 
  dplyr::rename(hv001 = dhsclust)
sf::st_geometry(coords) <- NULL


#.............
# K-means Clustering of DRC
#.............
keda <- data.frame(k = seq(2, 200, by = 1))
keda$kmeans <- map(keda$k, function(k){return(kmeans(x = coords[,2:3], centers = k))})
keda$wss <- map(keda$kmeans, "withinss")
keda$totalwss <- map(keda$wss, function(x){return(sum(x))})


keda.df <- keda %>% 
  dplyr::select(c("k", "totalwss")) %>% 
  tidyr::unnest()

kesteda.plotObj <- keda.df %>% 
  tibble::as.tibble(.) %>% 
  ggplot() +
  geom_line(aes(x=k, y=totalwss)) + 
  geom_point(aes(x=k, y=totalwss)) +
  geom_vline(xintercept = 15, color = "red", linetype = 2, alpha = 0.8) +
  theme_minimal() + 
  ylab("Total Within-Cluster Sum of Squares") +
  xlab("K")
plot(kesteda.plotObj)


# K means of 8 appears to be the infelction point
k <- kmeans(coords[,2:3], 14)
drcpart <- cbind.data.frame(coords, k = k$cluster) 

drcpart.plotObj <- drcpart %>% 
  ggplot() + 
  geom_sf(data = DRCprov) + 
  geom_point(aes(x=longnum, y=latnum, color = factor(k))) +
  theme(legend.position = "none")
plot(drcpart.plotObj)


drcpart <- drcpart %>% 
  dplyr::select(c("hv001", "k")) %>% 
  dplyr::rename(kmeansk = k)

saveRDS(drcpart, "data/derived_data/kmean_drcpartitioned.rds")


