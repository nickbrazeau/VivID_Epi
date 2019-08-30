# IMPORTS and dependencies
library(tidyverse)
library(sf)
source("R/00-functions_basic.R")
tol <- 1e-3

# Notes on Lonely PSUs
# http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="adjust")

#...........................................................
# Antimalarial Cluster Usage
#...........................................................
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
kr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDKR61FL.rds")

# liftover drug function
# per document, missing goes to NO
missingliftover <- function(x){
  x <- ifelse(x == 9 | is.na(x), 0, x) # 9 is missing
  return(x)
}


denom <- kr %>% 
  dplyr::filter(b8 < 5) %>% # less than 5 years -- dhs calls for b19<60, but no b19 variable in kr, ir, hr... this should be sufficient as it is less than 5 years
  dplyr::filter(haven::as_factor(b5) == "yes") %>% # currently alive
  dplyr::filter(haven::as_factor(h22) == "yes") %>% # had fever in last two weeks
  dplyr::select(c(paste0("ml13", letters[1:8]), "v001", "v005", "v023")) %>% 
  dplyr::select(-c("ml13g")) %>% # coded as NA in cd2013
  dplyr::mutate(v005 = v005/1e6)

# clean up
denom[, grepl("ml13", colnames(denom))] <- lapply(denom[, grepl("ml13", colnames(denom))], 
                                                  missingliftover) %>% dplyr::bind_cols(.)
# check
sapply(denom, summary)
# add any use in
denom$anyatm = as.numeric( apply(denom[,grepl("ml13", colnames(denom))],
                                 1, function(x){return(any( x == 1))}) )
# Note, some individuals took multiple drugs. OK because small percent 36/1560

kdsrv_fvr <- denom %>% srvyr:::as_survey_design(ids = v001, 
                                                strata = v023, 
                                                weights = v005)

kdsrv_fvr_clst <- kdsrv_fvr  %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(v001) %>% 
  dplyr::summarise(
    n = srvyr::survey_total(count),
    # fansidar_cont_clst = srvyr::survey_mean(ml13a),
    # chloroquine_cont_clst = srvyr::survey_mean(ml13b),
    # amodiaquine_cont_clst = srvyr::survey_mean(ml13c),
    # quinine_cont_clst = srvyr::survey_mean(ml13d),
    # act_cont_clst = srvyr::survey_mean(ml13e),
    # otherartm_cont_clst = srvyr::survey_mean(ml13f),
    # other_cont_clst = srvyr::survey_mean(ml13h),
    anyatm_cont_clst = srvyr::survey_mean(anyatm)) %>% 
  dplyr::rename(hv001 = v001)  %>% 
  dplyr::mutate(hv001 = as.numeric(hv001)) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))



#-----------------------------------------------------------
# Not all clusters had kids with fever in the last 2 weeks,
# Need to impute 8 missing clusters
#-----------------------------------------------------------
# find missing clusters
clst.all <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds") %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::filter(hv103 == 1) %>% 
  dplyr::select(c("hv001", "latnum", "longnum", "geometry")) %>% 
  dplyr::filter(!duplicated(.))


kdsrv_fvr_clst.imp <- dplyr::left_join(clst.all, kdsrv_fvr_clst, by = "hv001") 

sum(is.na(kdsrv_fvr_clst.imp$anyatm_cont_clst))

#.................
# Look to see if KNN is reasonable
#.................
mssngclust <- kdsrv_fvr_clst.imp %>% 
  dplyr::filter(is.na(n)) %>% 
  dplyr::select(-c("latnum", "longnum"))

notmssngclust <- kdsrv_fvr_clst.imp %>% 
  dplyr::filter(!is.na(n)) %>% 
  dplyr::select(-c("latnum", "longnum"))

ggplot() + 
  geom_sf(data=notmssngclust, color= "black") + 
  geom_sf(data=mssngclust, color= "red", size = 2, alpha  = 0.5) 

# going to have border issues but that's OK

find_mean_actuse_five_clusters <- function(missingcluster.sf, knownclusters.sf){
  # obviously this is not a general function to be exported
  dist <- raster::pointDistance(p1 = sf::as_Spatial(missingcluster.sf),
                               p2 = sf::as_Spatial(knownclusters.sf),
                               lonlat = T)
  
  # find 5 nearby clusters
  dist.sorted.5 <- sort(dist)[1:5]
  nrbyclstrs <- which(dist %in% dist.sorted.5)
  # sanity check
  if(length(nrbyclstrs) != 5){
    stop("nearest clusters don't match up")
  }
  nrbyclstrs <- knownclusters.sf$hv001[nrbyclstrs]
  
  ret <- knownclusters.sf %>% 
    dplyr::filter(hv001 %in% nrbyclstrs) %>% 
    dplyr::summarise(
      anyatm_cont_clst = mean(anyatm_cont_clst)
    )
  
  sf::st_geometry(ret) <- NULL # don't need spatial points anymore
  
  missingcluster.sf$anyatm_cont_clst <- ret$anyatm_cont_clst
  
  sf::st_geometry(missingcluster.sf) <- NULL # don't need spatial points anymore...
  return(missingcluster.sf)
}


# run it all
mssngclust.list <- mssngclust %>% 
  base::split(1:nrow(.))


mssngclust <- lapply(mssngclust.list, find_mean_actuse_five_clusters, knownclusters.sf = notmssngclust) %>% 
  dplyr::bind_rows(.)


# Combine to useful parts
sf::st_geometry(notmssngclust) <- NULL # don't need spatial points anymore :)

kdsrv_fvr_clst.imp <- dplyr::bind_rows(notmssngclust, mssngclust) %>% 
  dplyr::arrange(hv001)



# NOW SCALE
kdsrv_fvr_clst.imp <- kdsrv_fvr_clst.imp %>% 
  dplyr::mutate(
    # fansidar_cont_scale_clst = my.scale(logit(fansidar_cont_clst, tol = tol), center = T, scale = T),
    # chloroquine_cont_scale_clst = my.scale(logit(chloroquine_cont_clst, tol = tol), center = T, scale = T),
    # amodiaquine_cont_scale_clst = my.scale(logit(amodiaquine_cont_clst, tol = tol), center = T, scale = T),
    # quinine_cont_scale_clst = my.scale(logit(quinine_cont_clst, tol = tol), center = T, scale = T),
    # act_cont_scale_clst = my.scale(logit(act_cont_clst, tol = tol), center = T, scale = T),
    # otherartm_cont_scale_clst = my.scale(logit(otherartm_cont_clst, tol = tol), center = T, scale = T),
    # other_cont_scale_clst = my.scale(logit(other_cont_clst, tol = tol), center = T, scale = T),
    anyatm_cont_logit_clst = logit(anyatm_cont_clst, tol = tol),
    anyatm_cont_logit_scale_clst = my.scale(anyatm_cont_logit_clst, center = T, scale = T)
  )  %>% 
  dplyr::select(-c("n"))





saveRDS(file = "data/derived_data/vividepi_kids_act_use_imputed.rds", object = kdsrv_fvr_clst.imp)




