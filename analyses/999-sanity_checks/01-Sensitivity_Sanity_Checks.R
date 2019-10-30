#...............................................................................................
# Purpose of this script is to provide sanity checks on 
# various covariates
#...............................................................................................

library(tidyverse)
dt <- readRDS("data/derived_data/vividepi_recode.rds")
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")


#...............................................................................................
# General Sanity Notes
#...............................................................................................

# note, no hiv testing in 70; 271/318 lost to contamination


#...............................................................................................
# Look at Sample Collection Times
#...............................................................................................
xtabs(~dt$hvdate_dtdmy)
xtabs(~dt$hvyrmnth_dtmnth)
# Note, sample collection times not uniformly distributed around country (e.g.) 
# Kinshasha sampled Early
ggplot() +
  geom_sf(data = DRCprov) + 
  geom_point(data = dt, aes(x=longnum, y=latnum, color = hvyrmnth_dtmnth))


#................................................
# Timing Affects Precipitation?
# Timing Affects Temperature?
#................................................

# precipitation is a big issue
# justification for coarsing the precip to 
# study time-frame to avoid issues of positivit 
dt %>% 
  ggplot() +
  geom_boxplot(aes(x=hvyrmnth_dtmnth, y = precip_lag_cont_scale_clst))


# temperature not as much
dt %>% 
  ggplot() +
  geom_boxplot(aes(x=hvyrmnth_dtmnth, y = temp_lag_cont_scale_clst))



#................................................
# Temperature and Altitude Collinearity (again)
#................................................
dt %>% 
  ggplot() + 
  geom_point(aes(x= alt_dem_cont_scale_clst, y = temp_lag_cont_scale_clst)) + 
  stat_smooth(aes(x= alt_dem_cont_scale_clst, y = temp_lag_cont_scale_clst), method = "lm", se = FALSE)

cor(dt$alt_dem_cont_scale_clst[!is.na(dt$temp_lag_cont_scale_clst)], 
    dt$temp_lag_cont_scale_clst[!is.na(dt$temp_lag_cont_scale_clst)])


ggplot() +
  geom_sf(data = DRCprov) + 
  geom_point(data = dt, aes(x=longnum, y=latnum, color = temp_lag_cont_scale_clst, size = alt_dem_cont_scale_clst))




#................................................
# Urban Residence
#................................................

ggplot() +
  geom_sf(data = DRCprov) + 
  geom_point(data = dt, aes(x=longnum, y=latnum, color = urban_rura, shape = urban_rural_pca_fctb_clst))



