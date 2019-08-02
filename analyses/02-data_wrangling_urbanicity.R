# IMPORTS and dependencies
library(tidyverse)
library(sf)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3

# Notes on Lonely PSUs
# http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
options(survey.lonely.psu="adjust")

# spatial from the DHS -- these are cluster level vars
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")  %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum))
sf::st_geometry(dt) <- NULL
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6
  )


ge <- sf::st_as_sf(readRDS(file = "data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
ge <- ge %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(hv001 = dhsclust) %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

DRCprov <- sf::st_as_sf(readRDS("data/map_bases/gadm/gadm36_COD_1_sp.rds"))


#..............................
#### A note on urbanicity ####
#..............................

# Potential (significant) misclassification bias in the DHS DRC-II coding of 
# urban vs. rural as has been noted here https://journals.sagepub.com/doi/10.1177/0021909617698842
# and can be seen by comparing hv025/026 with population density, light density, build, etc.

#.............
# Urban from DHS
#.............
levels(factor(haven::as_factor(dt$hv025))) # no missing
dt$hv025_fctb <- haven::as_factor(dt$hv025)
dt$hv025_fctb = forcats::fct_relevel(dt$hv025_fctb, "rural")
# # check
xtabs(~hv025 + hv025_fctb, data = dt, addNA = T)


#.............
# cluster degree of "build"
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2014
summary(dt$built_population_2014)
hist(dt$built_population_2014)
sum(dt$built_population_2014 < 0.01)
median( dt$built_population_2014[dt$built_population_2014 < 0.05] )
hist( dt$built_population_2014[dt$built_population_2014 < 0.05] )
# DECISION: Will use a logit transformation to get back to the real-line (and scale)
# large number of 0s
dt$built_population_2014_scale <- my.scale(logit(dt$built_population_2014, tol = tol), center = T, scale = T) # use logit to transform to real line
hist(dt$built_population_2014_scale)
summary(dt$built_population_2014_scale); sd(dt$built_population_2014_scale) # despite skew, scale seems to work

#.............
# cluster night-time light density
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2015
summary(dt$nightlights_composite)
hist(dt$nightlights_composite)
hist( dt$nightlights_composite[dt$nightlights_composite < 0.05] )
hist( dt$nightlights_composite[dt$nightlights_composite > 0.05] )
# large number of 0s (again)
dt$nightlights_composite_scale <- my.scale(log(dt$nightlights_composite + tol), center = T, scale = T)
hist(dt$nightlights_composite_scale)
summary(dt$nightlights_composite_scale); sd(dt$nightlights_composite_scale) # despite skew, scale seems to work


#.............
# cluster worldpop population-density estimate
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2015
summary(dt$all_population_count_2015)
hist(dt$all_population_count_2015)
hist( dt$all_population_count_2015[dt$all_population_count_2015 < 5e4] )
hist( dt$all_population_count_2015[dt$all_population_count_2015 > 5e4] )
# no 0s here but a lot of small pops
dt$all_population_count_2015_scale <- my.scale(log(dt$all_population_count_2015 + tol), center = T, scale = T)
hist(dt$all_population_count_2015_scale)
summary(dt$all_population_count_2015_scale); sd(dt$all_population_count_2015_scale) # scale here seems to compensate


#.............
# Accessibility from Weiss
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2015
summary(dt$travel_times_2015)
hist(dt$travel_times_2015)
dt <- dt %>% 
  dplyr::mutate(travel_times_2015_scale = my.scale(log(travel_times_2015 + tol), center = T, scale = T))

hist(dt$travel_times_2015_scale) # many, many 0s -- these are urban centers/places near big towns
hist(dt$travel_times_2015_scale) # standardization doesn't look as good, but should capture urban v. rural well
summary(dt$travel_times_2015_scale); sd(dt$travel_times_2015_scale) 




#.............
# return scaled vars with linker
#.............

urbanmat <- dt[!duplicated(dt$hv001), 
               c("hv001", "hv025", 
                 "built_population_2014", "nightlights_composite",
                 "all_population_count_2015", "travel_times_2015",
                 "built_population_2014_scale", "nightlights_composite_scale",
                 "all_population_count_2015_scale", "travel_times_2015_scale")] 

urbanmat <- urbanmat %>% 
  dplyr::rename(
    built_population_2014_cont_clst = built_population_2014,
    nightlights_composite_cont_clst =  nightlights_composite,
    all_population_count_2015_cont_clst =  all_population_count_2015,
    travel_times_2015_cont_clst =  travel_times_2015,
    built_population_2014_cont_scale_clst =  built_population_2014_scale,
    nightlights_composite_cont_scale_clst =  nightlights_composite_scale,
    all_population_count_2015_cont_scale_clst =  all_population_count_2015_scale,
    travel_times_2015_cont_scale_clst =  travel_times_2015_scale
  )



#.............
# write out
#.............
saveRDS(urbanmat, file = "data/derived_data/vividepi_urban_recoded.rds")
