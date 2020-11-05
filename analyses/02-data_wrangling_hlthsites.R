#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle distance to public Health Sites
# in the DRC as a proxy to healthcare accessibility
# Using OSRM to calculate road duration/distances
#----------------------------------------------------------------------------------------------------
library(tidyverse)
#..................................
# import data
#..................................
# Get cluster locations
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")
# drop observations with missing geospatial data 
ge <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::select(c("hv001", "longnum", "latnum")) %>% 
  dplyr::filter(!duplicated(.)) %>% 
  dplyr::rename(dhsclust = hv001) # legacy

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
# of  mean duration to hospital 
#-----------------------------------------------------------------
durations <- osrm::osrmTable(src = ge.osrm,
                             dst = hlthsites.harvard.drc.osrm,
                             measure = "duration") # in minutes

#......................
# viz out
#......................
DRCprov <- readRDS("~/Documents/GitHub/VivID_Epi/data/map_bases/vivid_DRCprov.rds")
mindur <- apply(durations$durations, 1, min)
mindur_df <- tibble::tibble(dhsclust = ge.osrm$dhsclust,
                            mindur = mindur)
mindur_df <- ge %>% 
  dplyr::left_join(., mindur_df, by = "dhsclust")

mindur_df %>% 
  ggplot() + 
  #geom_sf(data = DRCprov, color = "#737373", fill = "#525252") +
  geom_point(aes(x = longnum, y = latnum, color = mindur)) + 
  scale_color_viridis_c("Duration in Minutes \n to Nearest Hospital")

mindur_df %>% 
  dplyr::mutate(twohrbin = ifelse(mindur >= 120, 1, 0)) %>% 
  ggplot() + 
  #geom_sf(data = DRCprov, color = "#737373", fill = "#525252") +
  geom_point(aes(x = longnum, y = latnum, color = twohrbin)) + 
  scale_color_viridis_c("Duration in Minutes \n to Nearest Hospital")

# look at urban rural
dt.nosf <- dt
sf::st_geometry(dt.nosf) <- NULL
urban_rura <- dt.nosf %>% 
  dplyr::select(c("hv001", "urban_rura")) %>% 
  dplyr::filter(!duplicated(.)) %>% 
  dplyr::rename(dhsclust = hv001)

mindur_df %>% 
  dplyr::mutate(twohrbin = ifelse(mindur >= 120, 1, 0)) %>% 
  dplyr::left_join(., urban_rura, by = "dhsclust") %>% 
  ggplot() + 
  #geom_sf(data = DRCprov, color = "#737373", fill = "#525252") +
  geom_point(aes(x = longnum, y = latnum, color = twohrbin, shape = urban_rura)) 


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
dir.create("data/derived_data/", recursive = TRUE)
saveRDS(object = hlthdist, file = "data/derived_data/hlthdist_out_minduration.rds")
saveRDS(object = hlthsites.harvard.drc, file = "data/derived_data/hlthsites_harvard_drc.rds")
