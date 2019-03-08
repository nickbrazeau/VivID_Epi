#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and clean the various covariates in
# the CD2013, weather, and other datasets
# 
# Notes: All variables that I will use will either contain a "_fctb/m" or "_cont" to show
#        that I have manipulated/investigated that variable.
#        Men/Women recode combinations (i.e. ha in one and hb in other for same covariate)
#        will be combined to be hab##
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")


#--------------------------------------------------------------
# Section 1:Pulling data-map file for all recodes
#-------------------------------------------------------------- 
# this code is now located under sandbox/01-liftover_attempts.R
# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ -- note this is version 7
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")



# drop observations with missing geospatial data 
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

#--------------------------------------------------------------
# Section 2: Looking at recodes, manual data wrangle
#-------------------------------------------------------------- 
# The purist in me would like to do this in a way that was automatic
# however, it is going to be a lot of work with lifting over all the recodes
# and the covariates are not all the same names or numbers in the various recodes
# which makes for a ton of corner cases, that will take more time to optimize
# than to do by hand... 

#..........................
# Exposure of Interests
#..........................
# SURVEY CHARACTERISTICS/WEIGHTS
# 1. hv005, cluster level weights
# 2. hv006, Month of interview
# 3. hv007, Year of interview
# 4. hv016, day of interview
#
#
#
# COINFECTIONS/BIOMARKER VARIABLES
# 1. pfldh coinfection
# 2. po18s coinfection
# 3. HIV Status (HIV03)
# 4. Currently pregnant? (HA54)
# 5. Hemoglobin Level adjust for altitude and smoking (HA56 & HB56)
#
# SOCIOECOLOGICAL VARIABLES
# 1. Type of Residence i.e. Rural or Urban (HV025)
# 2. Drinking Water (categorical: HV201)
# 4. Main floor material (categorical: HV213)
# 5. Main wall material (categorical: HV214)
# 6. Main roof material (categorical: HV215)
# 7. Shares toilet with other household (HV225)
# 8. Wealth Index (HV270)
# 9. Biological Sex (HV104)
# 10. Age (HV105)
# 11. Highest year of education completed (continous: HV108) 
# 12. Highest grade of education completed (categorical: HV107) 
# 13. Number of Children Under Age of 5 (continuous: HV014)
# 14. Number of Household Members (continuous: HV009)
# 15. Owns livestock, herds, or farm animals (hv246)

# MALARIA-INTERVENTIONS
# 1. Person slept under an LLIN net (HML20)/Net by Insecticide 
# 
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)


#..........................................................................................
#                                         SURVEY CHARACTERISTICS
#..........................................................................................
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6,
                hiv05_wi = hiv05/1e6
  )
#.............
# dates
#.............
dt <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste0(hv016, "/", hv006, "/", hv007)),
                hvyrmnth_dtmnth = factor(paste0(lubridate::year(hvdate_dtdmy), "-", lubridate::month(hvdate_dtdmy)),
                                       levels = c( "2013-8", "2013-9", "2013-11", "2013-12", "2014-1", "2014-2"))
                )


#.............
# households?? hard to say
#............
dt <- dt %>% 
  dplyr::mutate(hhid_lvl = factor( gsub("        | ", "", hhid)))



#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
#.............
# pfldh/po18s
#.............
dt <- dt %>% 
  dplyr::mutate(
    pfldhct_cont_ind = pfctmean,
    po18sct_cont_ind = poctcrrct,
    pv18sct_cont_ind = pvctcrrct,
    pfldh_fctb_ind = factor(pfldh, levels=c("0", "1"), labels=c("fal-", "fal+")),
    po18s_fctb_ind = factor(po18s, levels=c("0", "1"), labels=c("ov-", "ov+")),
    pv18s_fctb_ind = factor(pv18s, levels=c("0", "1"), labels=c("viv-", "viv+"))
  )

#.............
# HIV
#.............
# levels(factor(haven::as_factor(dt$hiv03)))
# hivlabels <- names( attributes(dt$hiv03)$labels )
# not the hiv03 covariate only has hiv negate or hiv positive, going to drop to just hiv+ or hiv-

dt <- dt %>% 
  dplyr::mutate(
    hiv03_fctb_ind = forcats::fct_drop(haven::as_factor(hiv03)),
    hiv03_fctb_ind = forcats::fct_relabel(hiv03_fctb_ind, ~ gsub(" ", "", .x, fixed = TRUE)), 
    hiv03_fctb_ind = forcats::fct_relabel(hiv03_fctb_ind, ~ gsub("positive", "+", .x)), 
    hiv03_fctb_ind = forcats::fct_relabel(hiv03_fctb_ind, ~ gsub("negative", "-", .x))
  ) 



#.............
# hemoglobin
#.............
levels(factor(haven::as_factor(dt$ha56)))
levels(factor(haven::as_factor(dt$hb56)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing 
# confirm no missing sex and then can use this variable to distinguish ha56 and hv56
xtabs(~haven::as_factor(dt$hv104), addNA = T)

dt <- dt %>% 
  dplyr::mutate(hab56_cont_ind = ifelse(haven::as_factor(hv104) == "female", ha56, hb56),
                hab56_cont_ind = ifelse(hab56_cont_ind %in% c("997", "999"), NA, hab56_cont_ind),
                hab56_cont_ind = as.numeric(hab56_cont_ind)/10,
                hab56_cont_ind_scale = scale(hab56_cont_ind, center = T, scale = T)) # this becomes (x-mu)/sd




#..........................................................................................
#                                  SOCIOECOLOGICAL VARIABLES
#...........................................................................................
#.............
# Urban
#.............
levels(factor(haven::as_factor(dt$hv025))) # no missing
dt$hv025_fctb_clst <- haven::as_factor(dt$hv025)
dt$hv025_fctb_clst = forcats::fct_relevel(dt$hv025_fctb_clst, "rural")
xtabs(~dt$hv025 + dt$hv025_fctb_clst)


#.............
# cluster type of town (urbanicity)
#.............
levels(factor(haven::as_factor(dt$hv026))) 
dt$hv026_fctm_clst <- haven::as_factor(dt$hv026)
xtabs(~dt$hv026 + dt$hv026_fctm_clst) # no missing so can drop
dt <- dt %>% dplyr::mutate(
  hv026_fctm_clst = forcats::fct_drop(haven::as_factor(dt$hv026_fctm_clst)),
  hv026_fctm_clst = forcats::fct_rev(forcats::fct_reorder(.f = hv026_fctm_clst, .x = hv026_fctm_clst, .fun = length)))

#.............
# cluster degree of "build"
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2014
summary(dt$built_population_2014)
hist(dt$built_population_2014)
dt$built_population_2014_logit_cont_clust <- logit(dt$built_population_2014, tol = 1e-3) # tolerance of 1e-3, transform back to real-line with logit 
hist(dt$built_population_2014_logit_cont_clust) # still not normal dist (as expected because of all the 0s) but better

dt <- dt %>% 
  dplyr::mutate(built_population_2014_logit_cont_clust_scale = scale(built_population_2014_logit_cont_clust, center = T, scale = T)) # this becomes (x-mu)/sd

#.............
# cluster night-time light density
#.............
# see explanation in the DHS GC manual 
# NOTE, this is from 2015
summary(dt$nightlights_composite)
hist(dt$nightlights_composite)
dt$nightlights_composite_log_cont_clust <- log(dt$nightlights_composite + 1e-3) # tolerance of 1e-3, not many, many 0s
hist(dt$nightlights_composite_log_cont_clust) # many, many 0s now become tolerance

dt <- dt %>% 
  dplyr::mutate(nightlights_composite_log_cont_clust_scale = scale(nightlights_composite_log_cont_clust, center = T, scale = T)) # this becomes (x-mu)/sd

#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f but still missing factor from haven
sum(is.na(dt$hv104))
dt <- dt %>% 
  dplyr::mutate(hv104_fctb_ind = haven::as_factor(hv104),
                hv104_fctb_ind = if_else(hv104_fctb_ind != "missing", hv104_fctb_ind, factor(NA)),
                hv104_fctb_ind =  forcats::fct_drop(hv104_fctb_ind), 
                hv104_fctb_ind = forcats::fct_rev(forcats::fct_reorder(.f = hv104_fctb_ind, .x = hv104_fctb_ind, .fun = length))
  ) # female to default (b/c 0 and larger SE)

#.............
# pregnant
#.............
levels(factor(haven::as_factor(dt$ha54)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing -- but note that no/don't know are collapsed 

dt <- dt %>% 
  dplyr::mutate(
    ha54_fctb_ind = haven::as_factor(dt$ha54),
    ha54_fctb_ind = if_else(ha54_fctb_ind != "missing", ha54_fctb_ind, factor(NA)),
    ha54_fctb_ind =  forcats::fct_drop(ha54_fctb_ind), 
    ha54_fctb_ind = forcats::fct_relevel(ha54_fctb_ind, "no/don't know")
  ) # no is obvious default based on n


#.............
# age
#.............
hist(dt$hv105)
dt <- dt %>% 
  dplyr::mutate(hv105_cont_ind = ifelse(hv105 %in% c("97", "98", "99"), NA, hv105),
                hv105_cont_ind_scale = scale(hv105_cont_ind, center = T, scale = T))



#.............
# main floor, wall, roof
#.............
# recode to rudimentary or non-rudimentary like they did here 
# https://pdfs.semanticscholar.org/e290/cf81bdb182696505952f37d1c910db86925a.pdf
# check to see how collinear wealth index and housing construction are since
# housing construction is included in the wealth index


# floor
floor <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_floor_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv213 = haven::as_factor(hv213), 
                hv213 =  forcats::fct_drop(hv213))
dt <- dt %>%
  left_join(x=., y=floor, by="hv213") %>% 
  dplyr::mutate(hv213_fctb_ind = factor(floortype),
                hv213_fctb_ind = forcats::fct_relevel(hv213_fctb_ind, "non-rudimentary")
  )
xtabs(~dt$hv213 + dt$hv213_fctb_ind, addNA = T)


# wall
wall <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_wall_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv214 = haven::as_factor(hv214), 
                hv214 =  forcats::fct_drop(hv214))
dt <- dt %>%
  left_join(x=., y=wall, by="hv214") %>% 
  dplyr::mutate(hv214_fctb_ind = factor(walltype),
                hv214_fctb_ind = forcats::fct_relevel(hv214_fctb_ind, "non-rudimentary")
  )

xtabs(~dt$hv214 + dt$hv214_fctb_ind, addNA = T)



# roof
roof <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_roof_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv215 = haven::as_factor(hv215), 
                hv215 =  forcats::fct_drop(hv215))
dt <- dt %>%
  left_join(x=., y=roof, by="hv215") %>% 
  dplyr::mutate(hv215_fctb_ind = factor(rooftype),
                hv215_fctb_ind = forcats::fct_relevel(hv215_fctb_ind, "non-rudimentary")
  )

xtabs(~dt$hv215 + dt$hv215_fctb_ind, addNA = T)


#.............
# wealth index 
#.............
levels(factor(haven::as_factor(dt$hv270))) 
sum(is.na(dt$hv270)) # no missing wealth
dt <- dt %>% 
  dplyr::mutate(hv270_fctm_ind = haven::as_factor(hv270),
                hv270_fctm_ind = forcats::fct_rev(forcats::fct_reorder(.f = hv270_fctm_ind, .x = hv270_fctm_ind, .fun = length)), # poorest is largest
                hv270_fctb_ind = factor(if_else(hv270_fctm_ind %in% c("poorest", "poorer", "middle"), 
                                     "poor_b", "rich_b"), levels = c("poor_b", "rich_b"))
                )



#.............
# education grade (categorical)
#.............
# levels(factor(haven::as_factor(dt$hv106)))
# dt <- dt %>% 
#   dplyr::mutate(hv106_fctm = haven::as_factor(hv106),
#                 hv106_fctm = if_else(!hv106_fctm %in% c("don't know", "missing"), hv106_fctm, factor(NA)),
#                 hv106_fctm =  forcats::fct_drop(hv106_fctm) # ordinal so don't relevel based on size
#   )
# 
# 
# #.............
# # education grade 
# #.............
# levels(factor(haven::as_factor(dt$hv106)))

#.............
# years of education (continuous)
#.............
hist(dt$hv108)
dt <- dt %>% 
  dplyr::mutate(hv108_cont_ind = ifelse(hv108 %in% c("97", "98", "99"), NA, hv108),
                hv108_cont_ind_scale = scale(hv108_cont_ind, center = T, scale = T))


#.............
# ethnicity
#.............
# question is whether you are congolese
levels(factor(haven::as_factor(dt$s113a))) 
levels(factor(haven::as_factor(dt$sm113a))) 

table(factor(haven::as_factor(dt$s113a))) 
table(factor(haven::as_factor(dt$sm113a))) 


dt <- dt %>% 
  dplyr::mutate(hv113a = ifelse(haven::as_factor(hv104) == "female", s113a, sm113a),
                hv113a_fctb_ind = factor(hv113a, levels = c(0,1), labels = c("no", "yes")))


# question of ethnicity, use liftover
levels(factor(haven::as_factor(dt$s114))) 
levels(factor(haven::as_factor(dt$sm114))) 

table(factor(haven::as_factor(dt$s114))) 
table(factor(haven::as_factor(dt$sm114))) 


dt <- dt %>% 
  dplyr::mutate(hv114a = ifelse(haven::as_factor(hv104) == "female", s114, sm114))

ethnic <-  readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_ethnic_liftover.csv")

dt <- dt %>%
  left_join(x=., y=ethnic, by="hv114a") %>% 
  dplyr::mutate(hv114a_fctm_ind = factor(ethnicgroup),
                hv114a_fctm_ind = forcats::fct_rev(forcats::fct_reorder(.f = hv114a_fctm_ind, .x = hv114a_fctm_ind, .fun = length)))

xtabs(~dt$hv114a + dt$hv114a_fctm_ind, addNA = T)



#------------------------------------------
# Owns livestock, herds, or farm animals
#------------------------------------------
summary(dt$hv246)
table(dt$hv246) # 9 is missing

dt <- dt %>% 
  dplyr::mutate(
    hv246_fctb_ind = haven::as_factor(dt$hv246),
    hv246_fctb_ind = if_else(hv246_fctb_ind != 9, hv246_fctb_ind, factor(NA)),
    hv246_fctb_ind =  forcats::fct_drop(hv246_fctb_ind), 
    hv246_fctb_ind = forcats::fct_relevel(hv246_fctb_ind, "no")
  ) 
xtabs(~ dt$hv246 + dt$hv246_fctb_ind)

#------------------------------------------
# Occupation
#------------------------------------------
# Going to dichotomize as suspected outdoor versus indoor
# Will code "others" as NA because cannot discern, only 2 obs

dt$hv717 <- ifelse(haven::as_factor(dt$hv104) == "female", dt$v717, dt$mv717)
table(dt$v717)
table(dt$mv717)

occupation <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_occupation_liftover.csv")
occupation$hv717 <- factor(occupation$hv717)

dt <- dt %>% 
  dplyr::mutate(hv717 = haven::as_factor(hv717), 
                hv717 =  forcats::fct_drop(hv717))
dt <- dt %>%
  left_join(x=., y=occupation, by="hv717") %>% 
  dplyr::mutate(hv717_fctb_ind = factor(jobconditions),
                hv717_fctb_ind = forcats::fct_relevel(hv717_fctb_ind, "indoors"))

xtabs(~ dt$hv717 + dt$hv717_fctb_ind)



#------------------------------------------
# children under 5 number
#------------------------------------------
summary(dt$hv014) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv014_cont_ind = hv014,
                hv014_cont_ind_scale = scale(hv014_cont_ind, center = T, scale = T))

#------------------------------------------
# total household members
#------------------------------------------
summary(dt$hv009) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv009_cont_ind = hv009,
                hv009_cont_ind_scale = scale(hv009_cont_ind, center = T, scale = T))



#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................
#.............
# LLIN for INDIVIDUAL 
#.............
xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml20), addNA = T)
# there are 49 people that slept under ITN" but not LLIN and a few NAs 

xtabs(~haven::as_factor(dt$hml19) + haven::as_factor(dt$hml20), addNA = T)
# there are 122 people that slept under "ever treated net" but not LLIN

xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml19), addNA = T)
# there are 69 people that slept under "ever treated net" but not a ITN

# BASED on this pattern, am just going to consider LLIN

dt <- dt %>% 
  dplyr::mutate(hml20_fctb_ind = haven::as_factor(hml20)) # no missing here surprisingly

#.............
# LLIN-type of Inseciticide for INDIVIDUAL 
# Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
#.............
# read insecticide liftover table
insctcd <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_insecticide_liftover.csv")

dt <- dt %>% 
  dplyr::mutate(hml7 = haven::as_factor(hml7)) %>% 
  left_join(x=., y=insctcd, by="hml7") %>% 
  dplyr::mutate(insctcd_fctm_ind = factor(ifelse(hml20_fctb_ind == "no", "none", insctcd)),
                insctcd_fctm_ind = forcats::fct_relevel(insctcd_fctm_ind, "none")
  )


# sanity checks
xtabs(~dt$insctcd_fctm_ind + dt$hml20_fctb_ind, addNA=T)
xtabs(~dt$insctcd_fctm_ind + dt$hml7, addNA=T)
xtabs(~dt$hml20_fctb_ind + dt$hml7, addNA=T)







#..........................................................................................
#                                PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
#..........................................................................................

#.............
# Cluster-Level Altitude
#.............

dt <- dt %>% 
  dplyr::mutate(alt_dem_cont_clust = ifelse(alt_dem == 9999, NA, alt_dem), # note no missing (likely dropped with missing gps)
                alt_dem_fctb_clust = factor(
                  ifelse(alt_dem_cont_clust > median(alt_dem_cont_clust), "high", "low"),
                  levels = c("low", "high")),
                alt_dem_clust_log = log(alt_dem_cont_clust),
                alt_dem_clust_log_scale = scale(alt_dem_clust_log, center = T, scale = T)
  )

#.............
# Temperature 
#.............
# dt <- dt %>% 
#   dplyr::mutate(mean_temperature_2015_cont = ifelse(mean_temperature_2015 < 0, NA, mean_temperature_2015)) # these are clearly errors but aren't labelled as such unless separate map file?
# 
# #.............
# # Rainfall
# #.............
# dt <- dt %>% 
#   dplyr::mutate(rainfall_2015_cont = ifelse(rainfall_2015 < 0, NA, rainfall_2015)) # these are clearly errors but aren't labelled as such unless separate map file?
#   # however these got knocked out by the missing cluster gps coord 
# 





#..........................................................................................
#                               CLUSTER LEVEL VARIABLES
#..........................................................................................
# aggregate these by cluster and then merge in with GE/hv001

#.............
# LLIN Cluster Usage
#.............


#.............
# Antimalarial Cluster Usage
#.............
# chloroquine




#..........................................................................................
#                               Final Write Out
#..........................................................................................
keep <- colnames(dt)[ grepl("_cont|_fctm|_fctb|_wi|_dt|_lvl", colnames(dt)) ]

# keep these variables for mapping
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
keep <- c("hv001", "hivrecode_barcode", "pfldh", "pfctmean", "po18s", "poctcrrct", "pv18s", "pvctcrrct", 
          "original_platemnum", "hv002", "hv021", "hv022", "hv023", "hvidx", "shnprovin",
          colnames(ge)[colnames(ge) != "dhsclust"], # dhsclust is same as hv001
          keep)


dtkp <- dt[, keep]

saveRDS(dtkp, file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")


