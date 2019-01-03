#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and clean the various covariates in
# the CD2013, weather, and other datasets
# 
# Notes: All variables that I will use will either contain a "_fct" or "_cont" to show
#        that I have manipulated/investigated that variable.
#        Men/Women recode combinations (i.e. ha in one and hb in other for same covariate)
#        will be combined to be hab##
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions.R")


#--------------------------------------------------------------
# Section 1:Pulling map file for all recodes
#-------------------------------------------------------------- 
# this code is now located under sandbox/01-liftover_attempts.R
# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
load("~/Documents/GitHub/VivID_Epi/data/vividepi_raw.rda")
dt <- merge_pr_plsmdm_gemtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge)

# drop observations with missing geospatial data
dt <- dt %>% 
  dplyr::filter(!is.na(dt$dhscc))



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
# 11. Highest year of education completed (continous: HV107) 
# 12. Number of Children Under Age of 5 (continuous: )
# 
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
  dplyr::mutate(hv005_cont = hv005/1e6,
                hv005 = hv005/1e6,
                hiv05_cont = hiv05/1e6,
                hiv05 = hiv05/1e6
  )
#.............
# dates
#.............
dt <- dt %>% 
  dplyr::mutate(hvdate_fct = lubridate::dmy(paste0(hv016, "/", hv006, "/", hv007)))





#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
#.............
# pfldh/po18s
#.............
dt <- dt %>% 
  dplyr::mutate(
    pfldh_fct = factor(pfldh, levels=c("0", "1"), labels=c("fal+", "fal-")),
    po18s_fct = factor(po18s, levels=c("0", "1"), labels=c("ov+", "ov-"))
  )

#.............
# HIV
#.............
# levels(factor(haven::as_factor(dt$hiv03)))
# hivlabels <- names( attributes(dt$hiv03)$labels )
# not the hiv03 covariate only has hiv negate or hiv positive so good to go there

dt <- dt %>% 
  dplyr::mutate(
    hiv03_fct = haven::as_factor(dt$hiv03)
  )

#.............
# pregnant
#.............
levels(factor(haven::as_factor(dt$ha54)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing 

dt <- dt %>% 
  dplyr::mutate(
    ha54_fct = haven::as_factor(dt$ha54)
  )


#.............
# hemoglobin
#.............
levels(factor(haven::as_factor(dt$ha56)))
levels(factor(haven::as_factor(dt$hb56)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing 
dt <- dt %>% 
  dplyr::mutate(hab56_cont = ifelse(!is.na(ha56), ha56, ifelse(!is.na(hb56), ha56, "error")),
                hab56_cont = ifelse(hab56_cont %in% c("997", "999"), NA, hab56_cont),
                hab56_cont = as.numeric(hab56_cont)/10)




#..........................................................................................
#                                  SOCIOECOLOGICAL VARIABLES
#...........................................................................................
#.............
# Urban
#.............
levels(factor(haven::as_factor(dt$ha54)))
dt$hv025 <- haven::as_factor(dt$hv025)


#.............
# drinking water and non-drink water
#.............
levels(factor(haven::as_factor(dt$hv201)))

dt <- dt %>% 
  dplyr::mutate(hv201_fct = haven::as_factor(hv201), # floor
                hv201_fct = ifelse(hv201_fct == "missing", NA, hv201_fct)
  )


#.............
# main floor, wall, roof
#.............
dt <- dt %>% 
  dplyr::mutate(hv213_fct = haven::as_factor(hv213), # floor
                hv213_fct = ifelse(hv213_fct == "missing", NA, hv213_fct), 
                
                hv214_fct = haven::as_factor(hv214), # wall
                hv214_fct = ifelse(hv214_fct == "missing", NA, hv214_fct), 
                
                hv215_fct = haven::as_factor(hv215), # roof
                hv215_fct = ifelse(hv215_fct == "missing", NA, hv215_fct)
                # no missing in roof?
  )



#.............
# shares toilet/toilet type
#.............
levels(factor(haven::as_factor(dt$hv205)))
levels(factor(haven::as_factor(dt$hv225)))

dt <- dt %>% 
  dplyr::mutate(hv205_fct = haven::as_factor(hv205), # toilet facility
                hv205_fct = ifelse(hv205_fct == "missing", NA, hv205_fct), 
                
                hv225_fct = haven::as_factor(hv225), # share toilet yes no
                hv225_fct = ifelse(hv225_fct == "missing", NA, hv225_fct)
  )



#.............
# wealth index 
#.............
levels(factor(haven::as_factor(dt$hv270))) 
sum(is.na(dt$hv270)) # no missing wealth
dt <- dt %>% 
  dplyr::mutate(hv270_fct = haven::as_factor(hv270))


#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f
sum(is.na(dt$hv104))
dt <- dt %>% 
  dplyr::mutate(hv104_fct = haven::as_factor(hv104))



#.............
# age
#.............
hist(dt$hv105)
dt <- dt %>% 
  dplyr::mutate(hv105_cont = ifelse(hv105 %in% c("97", "98", "99"), NA, hv105))

#.............
# education
#.............
levels(factor(haven::as_factor(dt$hv106)))
dt <- dt %>% 
  dplyr::mutate(hv106_fct = haven::as_factor(hv106),
                hv106_fct = ifelse(hv106_fct %in% c("don't know", "missing"), NA, hv106_fct)
  )




#------------------------------------------
# children under 5 number
#------------------------------------------
summary(dt$hv014) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv014_cont = hv014)


#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................
#.............
# LLIN
#.............
xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml20), addNA = T)
# there are 53 people that slept under ITN" but not LLIN

xtabs(~haven::as_factor(dt$hml19) + haven::as_factor(dt$hml20), addNA = T)
# there are 122 people that slept under "ever treated net" but not LLIN

xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml19), addNA = T)
# there are 69 people that slept under "ever treated net" but not a ITN

# BASED on this pattern, am just going to consider LLIN

dt <- dt %>% 
  dplyr::mutate(hml20_fct = haven::as_factor(hml20))

#.............
# LLIN-type of Inseciticide
# Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
#.............
# read insecticide liftover table
insctcd <- readr::read_csv("internal_datamap_files/pr_insecticide_liftover.csv")

dt <- dt %>% 
  dplyr::mutate(hml7 = haven::as_factor(hml7)) %>% 
  left_join(x=., y=insctcd, by="hml7") %>% 
  dplyr::mutate(insctcd_fct = factor(ifelse(hml20_fct == "no", "none", insctcd)))

# sanity checks
xtabs(~dt$insctcd_fct + dt$hml20_fct, addNA=T)
xtabs(~dt$insctcd_fct + dt$hml7, addNA=T)
xtabs(~dt$hml20_fct + dt$hml7, addNA=T)

#..........................................................................................
#                                PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
#..........................................................................................

#.............
# Cluster-Level Altitude
#.............
dt <- dt %>% 
  dplyr::mutate(hv040_cont = ifelse(hv040 == 9999, NA, hv040))

#.............
# CLIMATE
#.............


#.............
# Open Street Map
#.............
# R OSM
# http://osmar.r-forge.r-project.org/
# https://cran.r-project.org/web/packages/OpenStreetMap/index.html
# https://openmaptiles.com/downloads/dataset/satellite/africa/congo-democratic-republic/#2.9/-7.11/17.76














#..........................................................................................
#                               Final Write Out
#..........................................................................................
save(dt, DRCprov, ge, file = "data/vividepi_recode.rda")


