## .................................................................................
## Purpose: Make Satscan data
##
## Notes: 
## .................................................................................
library(tidyverse)
library(srvyr)
source("R/00-functions_basic.R")

#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
ge <- readRDS("data/raw_data/dhsdata/VivIDge.RDS")

# need weighed counts
data <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(clstn = srvyr::survey_total(count, vartype = c("se")),
                   pv18sn = srvyr::survey_total(pv18s, vartype = c("se"))) %>% 
  dplyr::mutate(clstn = round(clstn),
                pv18sn = round(pv18sn))

# bring in coords 
data <- dplyr::left_join(data, ge) %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2]) %>% 
  dplyr::select(c("hv001", "longnum", "latnum", "clstn", "pv18sn"))

# save out
dir.create("analyses/08-nonparam_riskfacts/data_for_satscan")
saveRDS(data, "analyses/08-nonparam_riskfacts/data_for_satscan/pvdata_for_satscan.RDS")