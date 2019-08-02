source("R/00-functions_Ensemble_Wrapper.R")
source("R/00-functions_iptw.R")
source("R/00-make_null_IPTW_distribs_brownian.R")
source("R/00-my_IPTW_performance_measure_energy.R")
set.seed(48, "L'Ecuyer")
library(tidyverse)
library(mlr)
library(rslurm)

dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL

plot(dt$hlth ~ dt$precip_lag_cont_clst)
hist(dt$annual_precipitation_2015)
mp <- dt %>% 
  dplyr::select(c("hv001", "temp_lag_cont_clst", "precip_lag_cont_clst", "hvyrmnth_dtmnth")) %>% 
  dplyr::filter(!duplicated(.))

mp.pv18s <- dtsrvy %>% 
  dplyr::select(c("hv001", "pv18s")) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(prev = srvyr::survey_mean(pv18s, vartype = c("se", "ci")))
                   
mp <- left_join(mp, mp.pv18s)



mp %>% ggplot() + 
  geom_boxplot(aes(x = hvyrmnth_dtmnth, y = precip_lag_cont_clst))


ge <- sf::st_as_sf(readRDS("data/raw_data/dhsdata/datasets/CDGE61FL.rds")) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(hv001 = dhsclust)

drcmp <- dplyr::left_join(mp, ge)

ggplot() +
  geom_sf(data=DRCprov) +
  geom_point(data = drcmp, aes(x=longnum, y=latnum, color = hvyrmnth_dtmnth))




