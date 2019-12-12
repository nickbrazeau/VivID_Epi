#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle distance to public Health Sites
# in the DRC as a proxy to healthcare accessibility
# Using OSRM to calculate road duration/distances
#----------------------------------------------------------------------------------------------------

#..................................
# import data
#..................................
# Get cluster locations
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# drop observations with missing geospatial data 
ge <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::select(c("hv001", "longnum", "latnum"))

hlthsites.harvard.drc <- readxl::read_excel("data/raw_data/harvard_dataverse/Ouma_Okiro_Snow_Africa_Hospitals_Data.xlsx") %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  magrittr::set_colnames(gsub(pattern = " ", "_", colnames(.))) %>% 
  dplyr::filter(country == "Democratic Republic of Congo") %>%
  sf::st_as_sf(coords = c("long", "lat"), 
               crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))


#..................................
# SPIN up Docker and start server
# then add these options
#..................................
# https://github.com/rCarto/osrm
remotes::install_github("rCarto/osrm")
library(osrm)
options(osrm.server = "http://0.0.0.0:5000/", osrm.profile = "driving")

ge.osrm <- ge %>% 
  dplyr::select("dhsclust") 
rownames(ge.osrm) <- ge.osrm$dhsclust


hlthsites.harvard.drc.osrm <- hlthsites.harvard.drc %>% 
  dplyr::mutate(id = paste0("hlth", seq(1:nrow(.))),
                id = factor(id)) %>%
  dplyr::select(id)

rownames(hlthsites.harvard.drc.osrm) <- hlthsites.harvard.drc.osrm$id


#-----------------------------------------------------------------
# Going to consider healthcare access as a function
# of how mean duration to hospital 
#-----------------------------------------------------------------
#' @param x sf object; query
#' @param y sf objest; target
#' @param long.distancematrix distance matrix, long format; distances between target and query
#' @param crs coordinate string
#' @param urbanwidth km for gbuffer
#' @param searchwidth km for gbuffer


osrm_distances_by_bcircle <- function(x, y, long.distancematrix, crs,
                                   searchwidth = 50){
  
  ret <- matrix(NA, nrow = nrow(x), ncol = 3)
  x <- sf::st_transform(x, crs)
  y <- sf::st_transform(y, crs)
  
  for(i in 1:nrow(x)){
    # get centroid
    centroid <- x[i,]
    centroid <- sf::as_Spatial(centroid) 
    # making bounding circle
    bcircle <- rgeos::gBuffer(centroid, width = searchwidth*1e3)
    
    # interset query
    hlthsites.catchment <- sf::st_intersection(y, sf::st_as_sf(bcircle))
    
    # ERROR CATCH clusters without catchment area
    # just find the min distance 
    if(nrow(hlthsites.catchment) == 0){
      
      durations <- osrm::osrmTable(src = centroid,
                                   dst = y,
                                   measure = "duration") # in minutes
      
      # fill in matrix
      hv001 <- x[i, "dhsclust"]
      sf::st_geometry(hv001) <- NULL 
      ret[i,1] <- unlist( hv001 ) # all this work to get a vector of 1
      ret[i, 2] <- min(durations$durations)
      ret[i, 3] <- NA
      
    } else {
      # get OSRM durations as planned
      durations <- osrm::osrmTable(src = centroid,
                                   dst = hlthsites.catchment,
                                   measure = "duration") # in minutes
      
      # fill in matrix
      hv001 <- x[i, "dhsclust"]
      sf::st_geometry(hv001) <- NULL 
      ret[i,1] <- unlist( hv001 ) # all this work to get a vector of 1
      ret[i, 2] <- mean(durations$durations)
      ret[i, 3] <- sd(durations$durations)
    }
    
  }
  
  colnames(ret) <- c("hv001", "mean_duration", "sd_duration")
  return(as.data.frame(ret))
  
}



hlthdist <- osrm_distances_by_bcircle(x = ge.osrm,
                                      y = hlthsites.harvard.drc.osrm,
                                      crs = "+proj=utm +zone=34 +datum=WGS84 +units=m",
                                      searchwidth = 100)


#..........
# Note, cluter 469 cannot be resolved with osrm
# do 5 knn again average
#..........
missclust <- data.frame(dhsclust = hlthdist$hv001[ which(is.na(hlthdist$mean_duration)) ] ) %>% 
  dplyr::left_join(., y=ge, by = "dhsclust") %>% 
  sf::st_as_sf(.)
knownclust <- data.frame(dhsclust = hlthdist$hv001[ which(! is.na(hlthdist$mean_duration)) ] ) %>% 
  dplyr::left_join(., y=ge, by = "dhsclust") %>% 
  sf::st_as_sf(.)


dist <- sf::st_distance(x = missclust,
                        y = knownclust,
                        which = "Great Circle")

# find 5 nearby clusters
dist.sorted.5 <- sort(dist)[1:5]
nrbyclstrs <- hlthdist$hv001[ which(dist %in% dist.sorted.5) ]

#.................
# sanity check
#.................
sanityplotdf <- ge %>% 
  dplyr::mutate(
    lvl = ifelse(dhsclust == 469, "miss", ifelse(
      dhsclust %in% nrbyclstrs, "nrby", "clst"
    ))
  )

sanitylabeldf <- sanityplotdf %>% 
  dplyr::filter(lvl == "nrby")

ggplot() + 
  geom_jitter(data = sanityplotdf, 
              aes(x=longnum, y = latnum, color = factor(lvl))) + 
  ggrepel::geom_label_repel(data = sanitylabeldf, 
                            aes(x=longnum, y = latnum, label = dhsclust))


# get their averages
hlthdist$mean_duration[hlthdist$hv001 == 469] <- mean(hlthdist$mean_duration[hlthdist$hv001 %in% nrbyclstrs])



#..........
# Now write out
#..........
saveRDS(object = hlthdist, file = "data/derived_data/hlthdist_out_minduration.rds")
saveRDS(object = hlthsites.harvard.drc, file = "data/derived_data/hlthsites_harvard_drc.rds")
