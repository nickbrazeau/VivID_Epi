#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and clean the various covariates for our analyses
# 
# Notes: All variables that I will use will either contain a "_fctb/m" or "_cont" to show
#        that I have manipulated/investigated that variable.
#        Men/Women recode combinations (i.e. ha in one and hb in other for same covariate)
#        will be combined to be hab##
#        
#        https://dhsprogram.com/pubs/pdf/FR300/FR300.pdf
#        
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
source("R/00-functions_basic.R")
tol <- 1e-3

#--------------------------------------------------------------
# Section 1:Pulling data-map file for all recodes
#-------------------------------------------------------------- 
# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ -- note this is version 7
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")


# drop observations with missing geospatial data 
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::filter(hv103 == 1) # subset to de-facto https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm

#--------------------------------------------------------------
# Section 2: Looking at recodes, manual data wrangle
#-------------------------------------------------------------- 

#..........................
# Exposure of Interests
#..........................
# SURVEY CHARACTERISTICS/WEIGHTS
# 1. hv005, cluster level weights
# 2. hv002, households


# COINFECTIONS/BIOMARKER VARIABLES
# 1. pfldh coinfection ; (personal)
# 3. HIV Status (HIV03) ; (personal)
# 4. Anemia Level (HA57 & HB57)
# 5. Duffy phenotype (wetlab result -- nearly all Duffyneg, so won't be formally considered in model)

# SOCIOECOLOGICAL VARIABLES
# 1. Biological Sex (HV104)
# 2. Age (HV105)
# 3. Main floor material (categorical: HV213)
# 3. Main wall material (categorical: HV214)
# 3. Main roof material (categorical: HV215)
# 3 -> 4. Building material (recode)
# 5. Wealth Index (recode; base wealth is hv270) 
# 6. Highest year of education completed (continous: HV108) 
# 7. Livestock ownership (categorical: HV246)
# 8. Occupation (categorical: "hv717")
# 9. Number of Household Members (continuous: HV009)
# 10. Number of Children Under Age of 5 (continuous: HV014)


# MALARIA-INTERVENTIONS
# 1. ITN Use
# 2. Cluster level antimalarial use
# 
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)
# 2. Temparature (gc recode)
# 3. Precipation (gc recode)
# 4. Urbanicity (my recode)



#########################################################################################################
##############                             SURVEY CHARACTERISTICS                          ##############
#########################################################################################################
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6
  )


#.............
# households
#............
summary(dt$hv002) # looks clean if we assume that households are numbered 1 - 34 in each cluster
hs <- dt %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = length(hv002), 
                   housemax = max(hv002))

summary(hs$housemax) # looks reasonable by cluster

dt <- dt %>% 
  dplyr::mutate(houseid = factor(paste0(hv001, "_", hv002)))


#.............
# dates
#.............
dt <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")))

# NOTE, some clusters have survey start and end dates that are in two months 
# (eg boundaries aren't clean/coinciding with a month. Naturally). Given
# grouping by month, need to assign a clusters "month" on the majority of days 
# that were spent surveying that clusters

# clusters without clean boundaries
clst_mnth_bounds <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
  dplyr::mutate(mnth = lubridate::month(hvdate_dtdmy)) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(moremnths = length(unique(mnth))) %>% 
  dplyr::filter(moremnths > 1)

clst_mnth_bounds.assign <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
  dplyr::filter(hv001 %in% clst_mnth_bounds$hv001) %>% 
  dplyr::mutate(hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
  dplyr::group_by(hv001, hvyrmnth_dtmnth) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n == max(n)) %>% 
  dplyr::select(-c("n"))

sf::st_geometry(clst_mnth_bounds.assign) <- NULL

dt <- dt %>% 
  dplyr::left_join(x=., y = clst_mnth_bounds.assign, by = "hv001") %>% 
  dplyr::mutate(hvyrmnth_dtmnth = ifelse(is.na(hvyrmnth_dtmnth),
                                         paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"),
                                         hvyrmnth_dtmnth))

dates <- readr::read_csv("internal_datamap_files/pr_date_liftover.csv")
dt <- dt %>% 
  dplyr::left_join(x=., y=dates, by = "hvyrmnth_dtmnth") %>% 
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))

xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)



#########################################################################################################
##############                          INDIVIDUAL LEVEL VARIABLES                         ##############
#########################################################################################################
#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
#.............
# pfldh/po18s
#.............
summary(dt$pfldh)
summary(dt$pv18s)

dt <- dt %>% 
  dplyr::mutate(
    pfldhct_cont = pfctmean,
    pfldhct_cont_log = log(pfctmean + tol),
    pv18sct_cont = pvctcrrct,
    pv18sct_cont_log = log(pvctcrrct + tol),
    po18sct_cont = poctcrrct,
    po18sct_cont_log = log(poctcrrct + tol),
    pfldh_fctb = factor(pfldh, levels=c("0", "1"), labels=c("falneg", "falpos")),
    pv18s_fctb = factor(pv18s, levels=c("0", "1"), labels=c("vivneg", "vivpos")),
    po18s_fctb = factor(po18s, levels=c("0", "1"), labels=c("ovneg", "ovpos"))
  )

dt[, colnames(dt)[grepl("pv18s|pfldh|po18s", colnames(dt))] ] %>% 
  sapply(., summary) # looks clean, all NAs are consistent except Pf but that has to do with double call strategy
                     # Remember, CT values >cutoff (species dependent) are coded as NA in raw data. Retained here
  
#.............
# HIV
#.............
levels(factor(haven::as_factor(dt$hiv03)))
hivlabels <- names( attributes(dt$hiv03)$labels ) # not the hiv03 covariate only has hiv negate or hiv positive, going to drop to just hiv+ or hiv-

dt <- dt %>% 
  dplyr::mutate(
    hiv03_fctb = forcats::fct_drop(haven::as_factor(hiv03)),
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub(" ", "", .x, fixed = TRUE)), 
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub("positive", "pos", .x)), 
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub("negative", "neg", .x)),
    hiv03_fctb = relevel(hiv03_fctb, "hivneg")
  ) 

xtabs(~hiv03_fctb + hiv03, data = dt, addNA = T)


#.............
# hemoglobin -- in ha56 (ha for females) and hb56 (hb for males)
# has been adjusted for altitude and smoking. 
# Note, even with the adjustments there are some extreme 
# readings (e.g. hemoglobins >20 and <5) that are clinically concerning...
#.............
levels(factor(haven::as_factor(dt$ha56))) # need to divide by 10; var is separate column for men and women
levels(factor(haven::as_factor(dt$hb56)))
# confirm no missing sex and then can use this variable to distinguish ha56 and hv56
xtabs(~haven::as_factor(dt$hv104), addNA = T)

dt <- dt %>% 
  dplyr::mutate(hab56_cont = ifelse(haven::as_factor(hv104) == "female", ha56, hb56),
                hab56_cont = ifelse(hab56_cont %in% c("997", "999"), NA, hab56_cont),
                hab56_cont = as.numeric(hab56_cont)/10,
                hab56_cont_scale = my.scale(hab56_cont))



# check hemoglobin recode for WOMEN
summary(dt$ha56)
nrow(dt) - length(dt$ha56[dt$hv104 == 2]) # missing tracks to male observations
summary(dt$ha56[dt$hv104 == 2 & dt$ha56 != 999 & dt$ha56 != 996 & dt$ha56 != 995])
sum(dt$ha56 %in% c(997, 999)) # 24 missing women hbs (23 missing, 1 inconsistent)

# check visually
dt %>% 
  dplyr::mutate(ha56 = ha56/10,
                ha56 = ifelse(ha56 %in% c(997, 999), NA, ha56)) %>% 
  ggplot() +
  geom_point(aes(x=hab56_cont, y=ha56)) + ylim(c(0,25)) + xlim(c(0,25))
# note, 7527 "rows" missing which corresponds to the 7503 males + 24 NAs coded in the dataset

# check hemoglobin recode for MEN
summary(dt$hb56)
summary(dt$hb56[dt$hv104 == 1 & dt$hb56 != 999 & dt$hb56 != 997])
sum(dt$hb56 %in% c(997, 999)) # 23 missing men hbs, 1 inconsistent

# check visually
dt %>% 
  dplyr::mutate(hb56 = hb56/10,
                hb56 = ifelse(hb56 %in% c(997, 999), NA, hb56)) %>% 
  ggplot() +
  geom_point(aes(x=hab56_cont, y=hb56)) + ylim(c(0,25)) + xlim(c(0,25))
# note, 8400 "rows" missing which corresponds to the 8376 females + 24 NAs coded in the dataset
hist(dt$hab56_cont)
xtabs(~dt$hab56_cont, addNA = T)



#.............
# Anemia
#.............
levels(factor(haven::as_factor(dt$ha57))) 
levels(factor(haven::as_factor(dt$hb57)))
# confirm no missing sex and then can use this variable to distinguish ha56 and hv56
xtabs(~haven::as_factor(dt$hv104), addNA = T)
# no missing anemia (that isn't dhs-recode dependent)
xtabs(~haven::as_factor(dt$ha57), addNA = T)
xtabs(~haven::as_factor(dt$hb57), addNA = T)

dt <- dt %>% 
  dplyr::mutate(hab57_fctb = ifelse(haven::as_factor(hv104) == "female", ha57, hb57),
                hab57_fctb = ifelse(hab57_fctb == 9, NA, hab57_fctb),
                hab57_fctb = ifelse(hab57_fctb == 4, "no", ifelse(
                  hab57_fctb %in% c(1:3), "yes", NA)),
                hab57_fctb = factor(hab57_fctb, levels = c("no", "yes")),
                hab57_fctb = relevel(hab57_fctb, "yes") # anemia is protective
  )

# check recode
xtabs(~dt$hab57_fctb + haven::as_factor(dt$ha57) + haven::as_factor(dt$hv104), addNA = T)  
xtabs(~dt$hab57_fctb + haven::as_factor(dt$hb57) + haven::as_factor(dt$hv104), addNA = T)  


# Recoding Anemia for with increased sensitivity/less strict cutoff
dt <- dt %>% 
  dplyr::mutate(lowhb_fctb = ifelse(haven::as_factor(hv104) == "female" & hab56_cont < 12.5, "yes", 
                                     ifelse(haven::as_factor(hv104) == "male" & hab56_cont < 13.5, "yes", 
                                            "no")),
                lowhb_fctb = factor(lowhb_fctb, levels = c("no", "yes"))
  )
                                     

xtabs(~dt$hab57_fctb + dt$lowhb_fctb, addNA = T)


#.............
# Duffy Phenotype
#.............
# wet lab. basically all individuals except for 3 -- no formal modeling needed (since so little variation)


#..........................................................................................
#                               DEMOGRAPHIC/BEHAVIORAL VARIABLES
#...........................................................................................
#.............
# sex
#.............
levels(factor(haven::as_factor(dt$hv104))) # no missing m/f but still missing factor from haven
sum(is.na(dt$hv104))
dt <- dt %>% 
  dplyr::mutate(hv104_fctb = haven::as_factor(hv104),
                hv104_fctb = forcats::fct_recode(hv104_fctb, NULL = "missing"),
                hv104_fctb =  forcats::fct_drop(hv104_fctb), 
                hv104_fctb = forcats::fct_rev(forcats::fct_reorder(.f = hv104_fctb, .x = hv104_fctb, .fun = length))
  ) # female to default (b/c 0 and larger SE)


#.............
# age
#.............
dt <- dt %>% 
  dplyr::mutate(hab1_cont = ifelse(haven::as_factor(hv104) == "female", ha1, hb1),
                hab1_cont = ifelse(hab1_cont %in% c("97", "98", "99"), NA, hab1_cont),
                hab1_cont_scale = my.scale(hab1_cont, center = T, scale = T)) 

summary(dt$hab1_cont) 
plot(dt$ha1, dt$hab1_cont)
plot(dt$hb1, dt$hab1_cont)

#.............
# main floor, wall, roof
#.............
# recode to rudimentary or non-rudimentary per previous manuscript (PMID: 28222094)
# then recode house to modern or traditional (final covar)
# https://pdfs.semanticscholar.org/e290/cf81bdb182696505952f37d1c910db86925a.pdf

# floor
floor <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_floor_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv213 = haven::as_factor(hv213), 
                hv213 =  forcats::fct_drop(hv213))
dt <- dt %>%
  left_join(x=., y=floor, by="hv213") %>% 
  dplyr::mutate( hv213_liftover = factor(floortype) )

xtabs(~dt$hv213 + dt$hv213_liftover, addNA = T)


# wall
wall <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_wall_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv214 = haven::as_factor(hv214), 
                hv214 =  forcats::fct_drop(hv214))
dt <- dt %>%
  left_join(x=., y=wall, by="hv214") %>% 
  dplyr::mutate( hv214_liftover = factor(walltype) )

xtabs(~dt$hv214 + dt$hv214_liftover, addNA = T)


# roof
roof <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_roof_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv215 = haven::as_factor(hv215), 
                hv215 =  forcats::fct_drop(hv215))
dt <- dt %>%
  left_join(x=., y=roof, by="hv215") %>% 
  dplyr::mutate( hv215_liftover = factor(rooftype) )

xtabs(~dt$hv215 + dt$hv215_liftover, addNA = T)

# Final Liftover to a binary variable
# 6 missing in floor; 5 missing in walls
wallfloormiss <- dt %>% 
  dplyr::select(c(hv213, hv214)) %>% 
  dplyr::mutate(hv213 = haven::as_factor(dt$hv213),
                hv214 = haven::as_factor(dt$hv214)) %>% 
  dplyr::filter(hv213 == "missing" | hv214 == "missing")
xtabs(~wallfloormiss$hv213 + wallfloormiss$hv214, addNA = T)
# different observations missing for floor v. wall


# Based on recent evidence, I think all metal roofs should be considered modern 
# because of the indoor temperature causing mosquito death PMC6533302
wallfloorroofrecode <- dt %>% 
  dplyr::select(c("hv213_liftover", "hv214_liftover", "hv215_liftover"))
sf::st_geometry(wallfloorroofrecode) <- NULL
dt <- dt %>% 
  dplyr::mutate(
    housecount = apply(wallfloorroofrecode, 1, function(x){return(sum(x == "non-rudimentary"))}),
    hv21345_fctb = ifelse(housecount == 3, "modern", "traditional"), # per PMID: 28222094
    hv21345_fctb = ifelse(haven::as_factor(dt$hv215) == "metal", "modern", hv21345_fctb), # per PMC6533302
    hv21345_fctb = factor(hv21345_fctb),
    hv21345_fctb = relevel(hv21345_fctb, "modern"))


# check -- seems reasonable
xtabs(~dt$hv21345_fctb, addNA = T)


#.............
# wealth index 
#.............
# As desrcibed in this manuscript (PMID: 28222094)
# housing materials are taken into consideration for wealth
# need to recode the wealth variable to avoid controlling for part of our
# effect when considering the covar housing materials (as they did above)
wlth <- readRDS(file = "data/derived_data/vividepi_wealth_recoded.rds")
dt <- dplyr::left_join(dt, wlth, by = "hivrecode_barcode")
xtabs(~dt$wlthrcde_fctm + haven::as_factor(dt$hv270)) # looks OK. Some poor/poorest got placed higher than expected
xtabs(~dt$wlthrcde_fctb + haven::as_factor(dt$hv270)) # looks OK. Some poor/poorest got placed higher than expected

dt <- dt %>% 
  dplyr::mutate(wlthrcde_fctb = relevel(wlthrcde_fctb, "not poor"))
  
#.............
# years of education (continuous)
#.............
haven::as_factor(dt$hv108)
dt <- dt %>% 
  dplyr::mutate(hv108_cont = ifelse(hv108 %in% c(97,98,99), NA, hv108),
                hv108_cont_scale = my.scale(hv108_cont))
#.............
# years of education (categorical)
#.............
edu <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_education_liftover.csv")
xtabs(~haven::as_factor(dt$hv106))
dt <- dt %>% 
  dplyr::mutate(hv106 = haven::as_factor(hv106), 
                hv106 =  forcats::fct_drop(hv106))
dt <- dt %>%
  left_join(x=., y=edu, by="hv106") %>% 
  dplyr::mutate( hv106_fctb = factor(hv106_fctb, levels = c("lower", "higher")),
                 hv106_fctb = relevel(hv106_fctb, "higher"))

xtabs(~dt$hv106 + dt$hv106_fctb, addNA = T)

#------------------------------------------
# Occupation as farmer or not
#------------------------------------------
table(factor(haven::as_factor(dt$hv104)), useNA = "always") # no missing m/f but still missing factor from haven
dt <- dt %>% 
  dplyr::mutate(occupation = ifelse(haven::as_factor(hv104) == "female",
                                    haven::as_factor(v717), 
                                    haven::as_factor(mv717)))

dt <- dt %>% 
  dplyr::mutate(farmer_fctb = ifelse(occupation %in% c(4,5), "farmer", "not farmer"),
                farmer_fctb = factor(farmer_fctb, levels = c("not farmer", "farmer"))) # not being a farmer protective
# note, we have coded missing as not a farmer

#------------------------------------------
# Owns livestock, herds, or farm animals
#------------------------------------------
summary(dt$hv246)
table(dt$hv246) # 9 is missing

dt <- dt %>% 
  dplyr::mutate(
    hv246_fctb = haven::as_factor(dt$hv246),
    hv246_fctb = forcats::fct_recode(hv246_fctb, NULL = "missing"),
    hv246_fctb =  forcats::fct_drop(hv246_fctb), 
    hv246_fctb = relevel(hv246_fctb, "no")
  ) 
xtabs(~ dt$hv246 + dt$hv246_fctb, addNA = T)

#------------------------------------------
# children under 5 number
#------------------------------------------
summary(dt$hv014) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv014_cont = hv014,
                hv014_cont_scale = my.scale(hv014_cont, center = T, scale = T))

#------------------------------------------
# total household members
#------------------------------------------
summary(dt$hv009) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv009_cont = hv009,
                hv009_cont_scale = my.scale(hv009_cont, center = T, scale = T))


#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................

#.............
# ITN for INDIVIDUAL 
#.............
# https://dhsprogram.com/Data/Guide-to-DHS-Statistics/index.cfm
# Use of Mosquito Nets by Persons in the Household
# Going to use the definition by Tusting et al. 2017 PMC5319641

xtabs(~haven::as_factor(dt$hml12)) # note, we have no one who slept under a both treated and untreated net
# and based on above, we have no missing net information

xtabs(~haven::as_factor(dt$hml12)) 
xtabs(~haven::as_factor(dt$hml10))

# using PMC5319641 definition
dt <- dt %>% 
  dplyr::mutate(
    ITN_fctb = ifelse(
      # (i) long-lasting insecticidal nets that were <= 3 y old at the time of survey
      !(haven::as_factor(hml4) %in% c("more than 3 years ago", "don't know", "missing")) & haven::as_factor(hml20) == "yes", 1,
      # (ii) conventional ITNs that were 1 y old 
      ifelse(haven::as_factor(hml4) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1, 
             # or were retreated within the year before the survey
             ifelse(haven::as_factor(hml9) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1, 
                    0))), # we know no missing net from above. they either reported yes or no at some level
    ITN_fctb = factor(ITN_fctb, levels = c(0,1), labels = c("no", "yes")),
    ITN_fctb = forcats::fct_drop(ITN_fctb),
    ITN_fctb = relevel(ITN_fctb, "yes")
  )


# sanity check
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml10)) 
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml12)) 
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml20)) 


#.............
# LLIN-Net
#.............
dt <- dt %>% 
  dplyr::mutate(hml20_fctb = haven::as_factor(hml20))


#.............
# LLIN-type of Inseciticide for INDIVIDUAL
# Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
#.............
# read insecticide liftover table
insctcd <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_insecticide_liftover.csv")

dt <- dt %>%
  dplyr::mutate(hml7 = haven::as_factor(hml7)) %>%
  left_join(x=., y=insctcd, by="hml7") %>%
  dplyr::mutate(insctcd_fctm = factor(ifelse(haven::as_factor(hml20) == "no", "none", insctcd)),
                insctcd_fctm = forcats::fct_relevel(insctcd_fctm, "none")
  )


# sanity checks
xtabs(~dt$insctcd_fctm + haven::as_factor(dt$hml20), addNA=T)
xtabs(~dt$insctcd_fctm + dt$hml7, addNA=T)
xtabs(~dt$hml20_fctb + dt$hml7, addNA=T)


# is this confounded by age of net?
xtabs(~dt$hml4 + dt$insctcd_fctm, addNA=T) # nope



#########################################################################################################
##############                             CLUSTER LEVEL VARIABLES                         ##############
#########################################################################################################
# ALL CLUSTER LEVEL VARIABLES WILL BE APPENDED WTIH "_clst"  
dtsrvy <- makecd2013survey(survey = dt)

#.............
# Weather
#.............
wthrnd.mth <- readRDS(file = "data/derived_data/vividep_weather_recoded_monthly.rds")
dt <- dt %>% 
  dplyr::left_join(., wthrnd.mth, by = c("hv001", "hvyrmnth_dtmnth_lag")) %>% 
  dplyr::mutate(
    precip_lag_cont_scale_clst = my.scale(precip_lag_cont_clst, center = T, scale = T),
    temp_lag_cont_scale_clst = my.scale(temp_lag_cont_clst, center = T, scale = T)
  )

wthrnd.mean <- readRDS(file = "data/derived_data/vividep_weather_recoded_mean.rds")
dt <- dt %>% 
  dplyr::left_join(., wthrnd.mean, by = c("hv001")) %>% 
  dplyr::mutate(
    precip_mean_cont_scale_clst = my.scale(precip_mean_cont_clst, center = T, scale = T),
    temp_mean_cont_scale_clst = my.scale(temp_mean_cont_clst, center = T, scale = T)
  )

#.............
# Cluster-Level Altitude
#.............
dt <- dt %>% 
  dplyr::mutate(alt_dem_cont_clst = ifelse(alt_dem == 9999, NA, alt_dem), # note no missing (likely dropped with missing gps)
                alt_dem_cont_scale_clst = my.scale(alt_dem_cont_clst, center = T, scale = T)
  )


#.............
# Distance to Water Source
#.............
wtrdist_out <- readRDS("data/derived_data/hotosm_waterways_dist.rds")
dt <- dt %>% 
  dplyr::left_join(x=., y = wtrdist_out, by = "hv001") %>% 
  dplyr::mutate(wtrdist_cont_log_clst = log(wtrdist_cont_clst + tol),
                wtrdist_cont_log_scale_clst = my.scale(wtrdist_cont_log_clst, center = T, scale = T)
                )


#..........................................................................................
#                           DEMOGRAPHIC/BEHAVIORAL VARIABLES
#..........................................................................................
#.............
# Urbanicity
#.............
dt <- dt %>% 
  dplyr::mutate(urban_rura_fctb = haven::as_factor(urban_rura))

#.............
# Health Care Access (Mean Distance to Health Site)
#.............
hlthdist_out <- readRDS("data/derived_data/hlthdist_out_minduration.rds") %>% 
  magrittr::set_colnames(c("hv001", "hlthdist_mean_duration", "hlthdist_sd_duration"))

dt <- dt %>% 
  dplyr::left_join(x=., y = hlthdist_out, by = "hv001") %>% 
  dplyr::mutate(
    hlthdist_cont_log_clst = log(hlthdist_mean_duration),
    hlthdist_cont_log_scale_clst = my.scale(hlthdist_cont_log_clst, center = T, scale = T),
    hlthst_duration_fctb_clst = ifelse(hlthdist_mean_duration > 60, "far", "near"),
    hlthst_duration_fctb_clst = factor(hlthst_duration_fctb_clst, levels = c("near", "far"))
    
  )




#.............
# Antimalarial Cluster Usage
#.............
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
actuse <- readRDS(file = "data/derived_data/vividepi_kids_act_use_imputed.rds")
dt <- dplyr::left_join(dt, actuse, by = "hv001")
dt %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    actuse = mean(anyatm_cont_clst),
    actuse.scale = mean(anyatm_cont_logit_scale_clst)) # looks good




#..........................................................................................
#                               Final Write Out
#..........................................................................................
saveRDS(dt, file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")


