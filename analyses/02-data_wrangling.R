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
# Section 1:Pulling data-map file for all recodes
#-------------------------------------------------------------- 
# this code is now located under sandbox/01-liftover_attempts.R
# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
load("~/Documents/GitHub/VivID_Epi/data/vividepi_raw.rda")
load("~/Documents/GitHub/VivID_Epi/data/vividspace_raw.rda")
dt <- merge_pr_plsmdm_ge_gc_mtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge, gc = gc)

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
  dplyr::mutate(hvdate_cont = lubridate::dmy(paste0(hv016, "/", hv006, "/", hv007)),
                hvyrmnth_fctm = factor(paste0(lubridate::year(hvdate_cont), "-", lubridate::month(hvdate_cont)),
                                       levels = c( "2013-8", "2013-9", "2013-11", "2013-12", "2014-1", "2014-2"))
                )


#.............
# households?? hard to say
#............
dt <- dt %>% 
  dplyr::mutate(hhid_fctm = factor( gsub("        | ", "", hhid)))



#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
#.............
# pfldh/po18s
#.............
dt <- dt %>% 
  dplyr::mutate(
    pfldh_fctb = factor(pfldh, levels=c("0", "1"), labels=c("fal-", "fal+")),
    po18s_fctb = factor(po18s, levels=c("0", "1"), labels=c("ov-", "ov+")),
    pv18s_fctb = factor(pv18s, levels=c("0", "1"), labels=c("viv-", "viv+"))
  )

#.............
# HIV
#.............
# levels(factor(haven::as_factor(dt$hiv03)))
# hivlabels <- names( attributes(dt$hiv03)$labels )
# not the hiv03 covariate only has hiv negate or hiv positive so good to go there, going to drop to just hiv+ or hiv-

dt <- dt %>% 
  dplyr::mutate(
    hiv03_fctb = forcats::fct_drop(haven::as_factor(dt$hiv03)),
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub(" ", "", .x, fixed = TRUE)), 
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub("positive", "+", .x)), 
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub("negative", "-", .x))
  ) 



#.............
# hemoglobin
#.............
levels(factor(haven::as_factor(dt$ha56)))
levels(factor(haven::as_factor(dt$hb56)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing 
dt <- dt %>% 
  dplyr::mutate(hab56_cont = ifelse(!is.na(ha56), ha56, ifelse(!is.na(hb56), ha56, "error")), # here I don't want ifelse to perserve type so I can check that I didn't make an error during import (i.e. can't coerce error to numeric)
                hab56_cont = ifelse(hab56_cont %in% c("997", "999"), NA, hab56_cont),
                hab56_cont = as.numeric(hab56_cont)/10)




#..........................................................................................
#                                  SOCIOECOLOGICAL VARIABLES
#...........................................................................................
#.............
# Urban
#.............
levels(factor(haven::as_factor(dt$hv025))) # no missing
dt$hv025_fctb <- haven::as_factor(dt$hv025)


#.............
# drinking water and non-drink water
#.............
levels(factor(haven::as_factor(dt$hv201)))

dt <- dt %>% 
  dplyr::mutate(hv201_fctm = haven::as_factor(hv201), # floor
                hv201_fctm = if_else(hv201_fctm != "missing", hv201_fctm, factor(NA)), # not sure why hv201_fctm == "missing", factor(NA), hv201_fctm misbehaves...
                hv201_fctm = forcats::fct_drop(hv201_fctm), # drop missing var we just did and any others that aren't supported
                hv201_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv201_fctm, .x = hv201_fctm, .fun = length)) # sort by that factors length, and then reverse
  )

# reproducible example of odd if_else behavior...
# set.seed(1)
# x <- factor(sample(letters[1:5], 10, replace = TRUE))
# x
# s <- if_else(x %in% c("a", "b", "c"), x, factor(NA))
# levels(s)
# xtabs(~s, addNA = T)
# 
# s <- if_else(x == "e", factor(NA), x)
# levels(s)
# xtabs(~s, addNA = T)



#.............
# main floor, wall, roof
#.............
dt <- dt %>% 
  dplyr::mutate(hv213_fctm = haven::as_factor(hv213), # floor
                hv213_fctm = if_else(hv213_fctm != "missing", hv213_fctm, factor(NA)), 
                hv213_fctm =  forcats::fct_drop(hv213_fctm), 
                hv213_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv213_fctm, .x = hv213_fctm, .fun = length)),
                
                hv214_fctm = haven::as_factor(hv214), # wall
                hv214_fctm = if_else(hv214_fctm != "missing", hv214_fctm, factor(NA)), 
                hv214_fctm =  forcats::fct_drop(hv214_fctm), 
                hv214_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv214_fctm, .x = hv214_fctm, .fun = length)),
                
                hv215_fctm = haven::as_factor(hv215), # roof
                hv215_fctm = if_else(hv215_fctm != "missing", hv215_fctm, factor(NA)), 
                hv215_fctm =  forcats::fct_drop(hv215_fctm), 
                hv215_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv215_fctm, .x = hv215_fctm, .fun = length))
                # no missing in roof?
  )



#.............
# shares toilet/toilet type
#.............
levels(factor(haven::as_factor(dt$hv205)))
levels(factor(haven::as_factor(dt$hv225)))

dt <- dt %>% 
  dplyr::mutate(hv205_fctm = haven::as_factor(hv205), # toilet facility
                hv205_fctm = if_else(hv205_fctm != "missing", hv205_fctm, factor(NA)), 
                hv205_fctm =  forcats::fct_drop(hv205_fctm), 
                hv205_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv205_fctm, .x = hv205_fctm, .fun = length)))
# dt <- dt %>%                
#  dplyr::mutate(hv225_fctb = haven::as_factor(hv225), # share toilet yes no
           #     hv225_fctb = if_else(hv225_fctb != "missing", hv225_fctb, factor(NA))
                # large number of people use bush/field (from hv205), so this variable is not independent. likely not adding much, dropping. 
#  )



#.............
# wealth index 
#.............
levels(factor(haven::as_factor(dt$hv270))) 
sum(is.na(dt$hv270)) # no missing wealth
dt <- dt %>% 
  dplyr::mutate(hv270_fctm = haven::as_factor(hv270),
                hv270_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv270_fctm, .x = hv270_fctm, .fun = length)) # poorest is largest
                )


#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f but still missing factor from haven
sum(is.na(dt$hv104))
dt <- dt %>% 
  dplyr::mutate(hv104_fctb = haven::as_factor(hv104),
                hv104_fctb = if_else(hv104_fctb != "missing", hv104_fctb, factor(NA)),
                hv104_fctb =  forcats::fct_drop(hv104_fctb), 
                hv104_fctb = forcats::fct_rev(forcats::fct_reorder(.f = hv104_fctb, .x = hv104_fctb, .fun = length))
                ) # female to default (b/c 0 and larger SE)

#.............
# pregnant
#.............
levels(factor(haven::as_factor(dt$ha54)))
# preglabels <- names( attributes(dt$ha54)$labels )
# not the ha54 covariate only has no/don't know or yes, so no extra missing -- but note that no/don't know are collapsed 

dt <- dt %>% 
  dplyr::mutate(
    ha54_fctb = haven::as_factor(dt$ha54),
    ha54_fctb = if_else(ha54_fctb != "missing", ha54_fctb, factor(NA)),
    ha54_fctb =  forcats::fct_drop(ha54_fctb), 
    ha54_fctb = forcats::fct_rev(forcats::fct_reorder(.f = ha54_fctb, .x = ha54_fctb, .fun = length))
  ) # no is obvious default based on n


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
  dplyr::mutate(hv106_fctm = haven::as_factor(hv106),
                hv106_fctm = if_else(!hv106_fctm %in% c("don't know", "missing"), hv106_fctm, factor(NA)),
                hv106_fctm =  forcats::fct_drop(hv106_fctm) # ordinal so don't relevel based on size
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
# there are 49 people that slept under ITN" but not LLIN and a few NAs 

xtabs(~haven::as_factor(dt$hml19) + haven::as_factor(dt$hml20), addNA = T)
# there are 122 people that slept under "ever treated net" but not LLIN

xtabs(~haven::as_factor(dt$hml10) + haven::as_factor(dt$hml19), addNA = T)
# there are 69 people that slept under "ever treated net" but not a ITN

# BASED on this pattern, am just going to consider LLIN

dt <- dt %>% 
  dplyr::mutate(hml20_fctb = haven::as_factor(hml20)) # no missing here surprisingly

#.............
# LLIN-type of Inseciticide
# Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
#.............
# read insecticide liftover table
insctcd <- readr::read_csv("internal_datamap_files/pr_insecticide_liftover.csv")

dt <- dt %>% 
  dplyr::mutate(hml7 = haven::as_factor(hml7)) %>% 
  left_join(x=., y=insctcd, by="hml7") %>% 
  dplyr::mutate(insctcd_fctm = factor(ifelse(hml20_fctb == "no", "none", insctcd)),
                insctcd_fctm = forcats::fct_relevel(insctcd_fctm, "none")
  )


# sanity checks
xtabs(~dt$insctcd_fctm + dt$hml20_fctb, addNA=T)
xtabs(~dt$insctcd_fctm + dt$hml7, addNA=T)
xtabs(~dt$hml20_fctb + dt$hml7, addNA=T)

#..........................................................................................
#                                PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
#..........................................................................................

#.............
# Cluster-Level Altitude
#.............
dt <- dt %>% 
  dplyr::mutate(hv040_cont = ifelse(hv040 == 9999, NA, hv040))


dt <- dt %>% 
  dplyr::mutate(mean_temperature_2015_cont = ifelse(mean_temperature_2015 < 0, NA, mean_temperature_2015)) # these are clearly errors but aren't labelled as such unless separate map file?

dt <- dt %>% 
  dplyr::mutate(rainfall_2015_cont = ifelse(rainfall_2015 < 0, NA, rainfall_2015)) # these are clearly errors but aren't labelled as such unless separate map file?
  # however these got knocked out by the missing cluster gps coord 






#..........................................................................................
#                               Final Write Out
#..........................................................................................
save(dt, file = "data/vividepi_recode.rda")


