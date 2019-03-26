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
devtools::install_github("malaria-atlas-project/malariaAtlas")
library(malariaAtlas)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3
set.seed(42)


bb <- osmdata::getbb("Democratic Republic of the Congo", 
                     featuretype = "country",
                     format_out = "sf_polygon")

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
# 2. hv006, Month of interview
# 3. hv007, Year of interview
# 4. hv016, day of interview
#
#
# COINFECTIONS/BIOMARKER VARIABLES
# 1. pfldh coinfection ; (personal + community)
# 2. po18s coinfection ; (personal + community)
# 3. HIV Status (HIV03) ; (personal + community)
# 4. Hemoglobin Level adjust for altitude and smoking (HA56 & HB56)
# 5. Duffy phenotype (wetlab work)

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
# 1. Person slept under an LLIN net (HML20) ; (personal + community)
# 2. Cluster level antimalarial use
# 
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)
# 2. Seasonality (ecological import)
# 3. Temparature (ecological import)
# 4. Precipation (ecological import)
# 5. Wetness Index (ecological import)
# 6. Type of Residence i.e. Rural or Urban (recode)
# 7. Ape Habitat Overlap (ecological import)



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
# dates
#.............
dt <- dt %>% 
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = factor(paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")),
                hvyrmnth_dtmnth_lag = factor(paste(lubridate::year(lubridate::rollback(hvdate_dtdmy, roll_to_first = F)), 
                                                   lubridate::month(lubridate::rollback(hvdate_dtdmy, roll_to_first = F)), sep = "-"))
  )
xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)

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
summary(dt$po18s)

dt <- dt %>% 
  dplyr::mutate(
    pfldhct_cont = pfctmean,
    pfldhct_cont_log = log(pfctmean + tol),
    pfldhct_cont_log_scale = scale(pfldhct_cont_log, center = T, scale = T),
    po18sct_cont = poctcrrct,
    po18sct_cont_log = log(poctcrrct + tol),
    po18sct_cont_log_scale = scale(po18sct_cont_log, center = T, scale = T),
    pv18sct_cont = pvctcrrct,
    pv18sct_cont_log = log(pvctcrrct + tol),
    pv18sct_cont_log_scale = scale(pv18sct_cont_log, center = T, scale = T),
    pfldh_fctb = factor(pfldh, levels=c("0", "1"), labels=c("falneg", "falpos")),
    po18s_fctb = factor(po18s, levels=c("0", "1"), labels=c("ovneg", "ovpos")),
    pv18s_fctb = factor(pv18s, levels=c("0", "1"), labels=c("vivneg", "vivpos")),
    pv18s_sens = ifelse(pv18sct_cont <= 42 & !is.na(pv18sct_cont), 1, 0), # empirical evidence to suggest 42 is an appropriate CT cutoff
    pv18s_fctb_sens = factor(pv18s_sens, levels=c("0", "1"), labels=c("vivneg", "vivpos")),
    pv18s_fctb_monoinfxn = factor(ifelse(pv18s_fctb == "vivpos" &  po18s_fctb == "ovneg" & pfldh_fctb == "falneg", 1, 0), 
                           levels=c("0", "1"), labels=c("vivmononeg", "vivmonopos"))
    
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
    hiv03_fctb = forcats::fct_relabel(hiv03_fctb, ~ gsub("negative", "neg", .x))
  ) 

xtabs(~hiv03_fctb + hiv03, data = dt, addNA = T)

#.............
# hemoglobin
#.............
levels(factor(haven::as_factor(dt$ha56))) # need to divide by 10; var is separate column for men and women
levels(factor(haven::as_factor(dt$hb56)))
# confirm no missing sex and then can use this variable to distinguish ha56 and hv56
xtabs(~haven::as_factor(dt$hv104), addNA = T)

dt <- dt %>% 
  dplyr::mutate(hab56_cont = ifelse(haven::as_factor(hv104) == "female", ha56, hb56),
                hab56_cont = ifelse(hab56_cont %in% c("997", "999"), NA, hab56_cont),
                hab56_cont = as.numeric(hab56_cont)/10,
                hab56_cont_scale = scale(hab56_cont, center = T, scale = T)) 

# check hemoglobin recode for WOMEN
    summary(dt$ha56)
    nrow(dt) - length(dt$ha56[dt$hv104 == 2]) # missing tracks
    summary(dt$ha56[dt$hv104 == 2 & dt$ha56 != 999 & dt$ha56 != 997])
    sum(dt$ha56 %in% c(997, 999)) # 25 missing women hbs
    # check visually
    dt %>% 
      dplyr::mutate(ha56 = ha56/10,
                    ha56 = ifelse(ha56 %in% c(997, 999), NA, ha56)) %>% 
      ggplot(.) +
      geom_point(aes(x=hab56_cont, y=ha56)) + ylim(c(0,25)) + xlim(c(0,25))
    # note, 7865 "rows" missing which corresponds to the 7840 males + 25 NAs coded in the dataset

# check hemoglobin recode for MEN
    summary(dt$hb56[dt$hv104 == 1 & dt$hb56 != 999 & dt$hb56 != 997])
    sum(dt$hb56 %in% c(997, 999)) # 25 missing women hbs
    # check visually
    dt %>% 
      dplyr::mutate(hb56 = hb56/10,
                    hb56 = ifelse(hb56 %in% c(997, 999), NA, hb56)) %>% 
      ggplot(.) +
      geom_point(aes(x=hab56_cont, y=hb56)) + ylim(c(0,25)) + xlim(c(0,25))
    # note, 8547 "rows" missing which corresponds to the 8523 females + 24 NAs coded in the dataset
hist(dt$hab56_cont_scale)
# final note on hemoglobin: a Hb > 18 is pretty high Likely due to 
# adjustment by alt and smoking (overadj).
# Few observations, going to assume negligible
# 49 missing. 
sum(dt$hab56_cont > 18, na.rm = T)
xtabs(~dt$hab56_cont, addNA = T)


#.............
# Duffy Phenotype
#.............
# TODO 


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
                hab1_cont = ifelse(hab56_cont %in% c("97", "98", "99"), NA, hab1_cont),
                hab1_cont_scale = scale(hab1_cont, center = T, scale = T)) 

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
floor <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_floor_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv213 = haven::as_factor(hv213), 
                hv213 =  forcats::fct_drop(hv213))
dt <- dt %>%
  left_join(x=., y=floor, by="hv213") %>% 
  dplyr::mutate( hv213_liftover = factor(floortype) )

xtabs(~dt$hv213 + dt$hv213_liftover, addNA = T)


# wall
wall <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_wall_liftover.csv")
dt <- dt %>% 
  dplyr::mutate(hv214 = haven::as_factor(hv214), 
                hv214 =  forcats::fct_drop(hv214))
dt <- dt %>%
  left_join(x=., y=wall, by="hv214") %>% 
  dplyr::mutate( hv214_liftover = factor(walltype) )

xtabs(~dt$hv214 + dt$hv214_liftover, addNA = T)


# roof
roof <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_roof_liftover.csv")
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

dt <- dt %>% 
  dplyr::mutate(hv21345_fctb = ifelse(apply(dt[,c("hv213_liftover", "hv214_liftover", "hv215_liftover")], 
                                            1, 
                                            function(x){
                                              return( all(x == "non-rudimentary", na.rm = F) ) # na.rm F to perserve missing
                                              }), 
                                      "modern", "traditional"),
                hv21345_fctb = factor(hv21345_fctb))
                
# check -- seems reasonabl
xtabs(~dt$hv21345_fctb)


#.............
# wealth index 
#.............
# As desrcibed in this manuscript (PMID: 28222094)
# housing materials are  taken into consideration for wealth
# need to recode the wealth variable to avoid controlling for part of our
# effect when considering the covar housing materials (as they did above)
wlth <- readRDS(file = "data/derived_data/vividepi_wealth_recoded.rds")
dt <- dplyr::left_join(dt, wlth, by = "hivrecode_barcode")
xtabs(~dt$wlthrcde_fctm + haven::as_factor(dt$hv270)) # looks OK. Some poor/poorest got placed higher than expected



#.............
# years of education (continuous)
#.............
hist(dt$hv108)
summary(dt$hv108)
dt <- dt %>% 
  dplyr::mutate(hv108_cont = ifelse(hv108 %in% c("97", "98", "99"), NA, hv108),
                hv108_cont_scale = scale(hv108_cont, center = T, scale = T))
xtabs(~dt$hv108 + dt$hv108_cont, addNA = T)


# #.............
# # ethnicity
# #.............
# A note on ethnicity: resolution not great here. Really care about Duffy status. Captured with
# genotyping. 
# 
# 
# # question is whether you are congolese
# levels(factor(haven::as_factor(dt$s113a))) 
# levels(factor(haven::as_factor(dt$sm113a))) 
# 
# table(factor(haven::as_factor(dt$s113a))) 
# table(factor(haven::as_factor(dt$sm113a))) 
# 
# 
# dt <- dt %>% 
#   dplyr::mutate(hv113a = ifelse(haven::as_factor(hv104) == "female", s113a, sm113a),
#                 hv113a_fctb = factor(hv113a, levels = c(0,1), labels = c("no", "yes")))
# 
# 
# # question of ethnicity, use liftover
# levels(factor(haven::as_factor(dt$s114))) 
# levels(factor(haven::as_factor(dt$sm114))) 
# 
# table(factor(haven::as_factor(dt$s114))) 
# table(factor(haven::as_factor(dt$sm114))) 
# 
# 
# dt <- dt %>% 
#   dplyr::mutate(hv114a = ifelse(haven::as_factor(hv104) == "female", s114, sm114))
# 
# ethnic <-  readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_ethnic_liftover.csv")
# 
# dt <- dt %>%
#   left_join(x=., y=ethnic, by="hv114a") %>% 
#   dplyr::mutate(hv114a_fctm = factor(ethnicgroup),
#                 hv114a_fctm = forcats::fct_rev(forcats::fct_reorder(.f = hv114a_fctm, .x = hv114a_fctm, .fun = length)))
# 
# xtabs(~dt$hv114a + dt$hv114a_fctm, addNA = T)
# 


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
    hv246_fctb = forcats::fct_relevel(hv246_fctb, "no")
  ) 
xtabs(~ dt$hv246 + dt$hv246_fctb, addNA = T)

#------------------------------------------
# Occupation
#------------------------------------------
# Going to dichotomize as suspected outdoor versus indoor
# Will code "others" as NA because cannot discern, only 2 obs

dt$hab717 <- ifelse(haven::as_factor(dt$hv104) == "female", dt$v717, dt$mv717)
table(dt$v717)
table(dt$mv717)
table(dt$hab717)

occupation <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_occupation_liftover.csv")
occupation$hab717 <- factor(occupation$hv717) # my mac being a pain

dt <- dt %>% 
  dplyr::mutate(hab717 = haven::as_factor(hab717), 
                hab717 =  forcats::fct_drop(hab717))
dt <- dt %>%
  left_join(x=., y=occupation, by="hab717") %>% 
  dplyr::mutate(hab717_fctb = factor(jobconditions),
                hab717_fctb = forcats::fct_relevel(hab717_fctb, "indoors"))

xtabs(~dt$hab717 + dt$hab717_fctb, addNA = T)

#------------------------------------------
# children under 5 number
#------------------------------------------
summary(dt$hv014) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv014_cont = hv014,
                hv014_cont_scale = scale(hv014_cont, center = T, scale = T))

#------------------------------------------
# total household members
#------------------------------------------
summary(dt$hv009) # looks clean
dt <- dt %>% 
  dplyr::mutate(hv009_cont = hv009,
                hv009_cont_scale = scale(hv009_cont, center = T, scale = T))


#..........................................................................................
#                                 MALARIA-INTERVENTIONS
#..........................................................................................
#.............
# Health Insurance
#.............
dt <- dt %>% 
  dplyr::mutate(hab481 =  ifelse(haven::as_factor(hv104) == "female", v481, mv481),
                hab481 = ifelse(hab481 == 9, NA, hab481))
xtabs(~dt$hab481, addNA = T)

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
summary(dt$hml20) # no missing here

dt <- dt %>% 
  dplyr::mutate(hml20_fctb = haven::as_factor(hml20))

#.............
# LLIN-type of Inseciticide for INDIVIDUAL
# Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
#.............
# read insecticide liftover table
insctcd <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_insecticide_liftover.csv")

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





#########################################################################################################
##############                             CLUSTER LEVEL VARIABLES                         ##############
#########################################################################################################
# ALL CLUSTER LEVEL VARIABLES WILL BE APPENDED WTIH "_clst"  
dtsrvy <- makecd2013survey(survey = dt)

#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
Hbmiss <- dt[is.na(dt$hab56_cont),]
table(Hbmiss$hv001) # looks well dispersed among clusters

coinfxnbiomrk <- dtsrvy %>% 
  dplyr::group_by(hv001) %>% # cluster level 
  dplyr::summarise(
    pfldh_cont_clst = srvyr::survey_mean(x = pfldh),
    po18s_cont_clst = srvyr::survey_mean(x = po18s),
    hiv03_cont_clst = srvyr::survey_mean(x = hiv03),
    hab56_cont_clst = srvyr::survey_quantile(hab56_cont, quantiles = c(0.5), vartype = c("se"), na.rm = T)) %>% 
  dplyr::mutate(
    pfldh_cont_scale_clst = scale(logit(pfldh_cont_clst, tol = tol), center = T, scale = T),
    po18s_cont_scale_clst = scale(logit(po18s_cont_clst, tol = tol), center = T, scale = T),
    hiv03_cont_scale_clst = scale(logit(hiv03_cont_clst, tol = tol), center = T, scale = T),
    hab56_cont_scale_clst = scale(hab56_cont_clst_q50, center = T, scale = T)
  ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))

sapply(coinfxnbiomrk, summary) # looks clean, had to remove NAs in Hb because 49 missing (spread across 39 clusters)

dt <- dplyr::left_join(x = dt, y = coinfxnbiomrk)


#..........................................................................................
#                                ECOLOGICAL VARIABLES
#..........................................................................................
source("R/00-functions_ecology.R")
#.............
# Seasonality
#.............




#.............
# Precipitation (CHRIPS) & Temperature (Emch/Manny)
#.............
# precip data
precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T)
precipfrst <- lapply(precip, readRasterBB, bb = bb)
precipdf <- tibble::tibble(names = basename(precip)) %>% 
  dplyr::mutate(names = gsub("chirps-v2.0.|.tif", "", names),
                year = stringr::str_split_fixed(names, "\\.", n=2)[,1] ,
                mnth =  stringr::str_split_fixed(names, "\\.", n=2)[,2] ,
                hvdate_dtdmy = lubridate::dmy(paste(1, mnth, year, sep = "/")),
                year = lubridate::year(hvdate_dtdmy),
                mnth = lubridate::month(hvdate_dtdmy),
                hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                precipraster = precipfrst) %>% 
  dplyr::select(c("hvyrmnth_dtmnth_lag", "precipraster"))

# temp data
temp <- raster::stack("data/raw_data/weather_data/emch_manny/cru_ts4.01.2011.2016.tmp.dat.nc") %>% 
        raster::crop(x = ., y = sf::as_Spatial(bb))

tempdf <- tibble::tibble(orignnames = names(temp),
                         names = gsub("X", "", names(temp)),
                         hvdate_dtdmy = lubridate::ymd(names),
                         year = lubridate::year(hvdate_dtdmy),
                         mnth = lubridate::month(hvdate_dtdmy),
                         hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
                         tempraster = lapply(orignnames, raster::subset, x = temp)) %>% 
  dplyr::select(c("hvyrmnth_dtmnth_lag", "tempraster"))



wthrnd <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>% 
  dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
wthrnd <- wthrnd[!duplicated(wthrnd),]

wthrnd <- wthrnd %>% 
  dplyr::left_join(., tempdf) %>% 
  dplyr::left_join(., precipdf)

# Drop in a for loop
wthrnd$precip_lag_cont_clst <- NA
wthrnd$temp_lag_cont_clst <- NA

for(i in 1:nrow(wthrnd)){
  # precip
  wthrnd$precip_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd$precipraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
                    )
  
  # temp
  wthrnd$temp_lag_cont_clst[i] <- 
    raster::extract(x = wthrnd$tempraster[[i]],
                    y = sf::as_Spatial(wthrnd$geometry[i]),
                    buffer = wthrnd$buffer[i],
                    fun = mean,
                    sp = F
    )
  
}

wthrnd <- wthrnd %>% 
  dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst", "temp_lag_cont_clst")) %>% 
  dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))

dt <- dt %>% 
  dplyr::left_join(., wthrnd, by = c("hv001", "hvyrmnth_dtmnth_lag")) %>% 
  dplyr::mutate(precip_lag_cont_log_clst = log(precip_lag_cont_clst + tol),
                temp_lag_cont_log_clst = log(temp_lag_cont_clst + tol))



#.............
# Cluster-Level Altitude
#.............
dt <- dt %>% 
  dplyr::mutate(alt_dem_cont_clst = ifelse(alt_dem == 9999, NA, alt_dem), # note no missing (likely dropped with missing gps)
                alt_dem_fctb_clst = factor(
                  ifelse(alt_dem_cont_clst > median(alt_dem_cont_clst), "high", "low"),
                  levels = c("low", "high")),
                alt_dem_clst_log = log(alt_dem_cont_clst + tol),
                alt_dem_clst_log_scale = scale(alt_dem_clst_log, center = T, scale = T)
  )

#.............
# Ape Habitat Overlap
#.............
drc_ape <- readRDS(file = "data/redlist_species_data_primate/drc_ape.rds")
drc_ape_simplified <- rgeos::gSimplify(
  sf::as_Spatial(sf::st_union(drc_ape)), 
  tol = 1e-3)

dt <- dt %>% 
  dplyr::mutate(ape_habitat_fctb_clst = factor(
    ifelse(rgeos::gContains(spgeom1 = drc_ape_simplified,
                            spgeom2 = sf::as_Spatial(sf::st_as_sf(.)),
                            byid = T),
           1, 0),
    levels = c(0,1), labels = c("NoOverlap", "Overlap"))
  )
xtabs(~dt$ape_habitat_fctb_clst)

#.............
# Distance to Water Source
#.............
# TODO 
# Use IDW for raster distance
# clsts <- dt[!duplicated(dt$hv001), c("hv001", "geometry", "urban_rura")]
# wtrdist <- geosphere::dist2Line(p = clsts, line = wtr, distfun = distGeo)


#..........................................................................................
#                           DEMOGRAPHIC/BEHAVIORAL VARIABLES
#..........................................................................................
#.............
# Urbanicity
#.............
#..............................
#### A note on urbanicity ####
#..............................
# Potential (significant) misclassification bias in the DHS DRC-II coding of 
# urban vs. rural as has been noted here https://journals.sagepub.com/doi/10.1177/0021909617698842
# and can be seen by comparing hv025/026 with population density, light density, build, and accessibility 
# the first PCA explains ~80% of the variation. Will use that as my new "urban score" 

urb <- readRDS(file = "data/derived_data/vividepi_urban_recoded.rds")
dt <- dplyr::left_join(dt, urb, by = "hv001")
xtabs(~dt$urbanscore_fctm_clust + haven::as_factor(dt$hv025)) # looks OK. Some rural places are now urban, which is consistent with what I expected

#.............
# Distance to Health Site
#.............
# TODO 
# htlhdist <- geosphere::distm(x = clsts, y=hlthoffc, fun = distGeo)


#.............
# LLIN Cluster Usage; median wealth; median educ; health insur
#.............
summary(dt$hml20) # no missing
summary(dt$wlthrcde_fctm)
summary(dt$hv108_cont)

democlust <- dtsrvy %>% 
  dplyr::mutate(wlthrcde_fctm_ord = factor(wlthrcde_fctm,
                                           levels = c("poorest", "poor", "middle", "rich", "richest"),
                                           ordered = T),
                wlthrcde_fctm_ord_num = as.numeric(wlthrcde_fctm_ord)) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    hml20_cont_clst = srvyr::survey_mean(hml20, vartype = c("se")),
    wlthrcde_fctm_clst = srvyr::survey_median(wlthrcde_fctm_ord_num, quantiles = c(0.5), vartype = c("se")),
    hv108_cont_clst = srvyr::survey_median(hv108_cont, quantiles = c(0.5), vartype = c("se"), na.rm = T),
    hab481_cont_clst = srvyr::survey_mean(hab481, vartype = c("se"), na.rm = T)
  ) %>% 
  dplyr::mutate(
    hml20_cont_scale_clst = scale(logit(hml20_cont_clst, tol = tol), center = T, scale = T),
    hv108_cont_scale_clst = scale(log(hv108_cont_clst_q50 + tol), center = T, scale = T),
    wlthrcde_fctm_clst_q50_fctm = factor(floor(wlthrcde_fctm_clst_q50), # taking lower bracket of wealth if median was split
                                         levels = c(1,2,3,4,5),
                                         labels = c("poorest", "poor", "middle", "rich", "richest")),
    hab481_cont_scale_clst = scale(logit(hab481_cont_clst, tol = tol), center = T, scale = T)
          )


sapply(democlust, summary) # looks clean

dt <- dplyr::left_join(x = dt, y = democlust)


#.............
# Antimalarial Cluster Usage
#.............
#https://dhsprogram.com/data/Guide-to-DHS-Statistics/
kr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDKR61FL.rds")
# liftover drug function
# per document, missing goes to NO
missingliftover <- function(x){
  x <- ifelse(x == 9 | is.na(x), 0, x) # 9 is missing
  return(x)
}

denom <- kr %>% 
  dplyr::filter(v012 < 60) %>% # less than 5 years
  dplyr::filter(haven::as_factor(b5) == "yes") %>% # currently alive
  dplyr::filter(haven::as_factor(h22) == "yes") %>% # had fever in last two weeks
  dplyr::select(c(paste0("ml13", letters[1:8]), "v001", "v005", "v023")) %>% 
  dplyr::select(-c("ml13g")) %>% 
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
    fansidar_cont_clst = srvyr::survey_mean(ml13a),
    chloroquine_cont_clst = srvyr::survey_mean(ml13b),
    amodiaquine_cont_clst = srvyr::survey_mean(ml13c),
    quinine_cont_clst = srvyr::survey_mean(ml13d),
    act_cont_clst = srvyr::survey_mean(ml13e),
    otherartm_cont_clst = srvyr::survey_mean(ml13f),
    other_cont_clst = srvyr::survey_mean(ml13h),
    anyatm_cont_clst = srvyr::survey_mean(anyatm)) %>% 
  dplyr::mutate(
    fansidar_cont_scale_clst = scale(logit(fansidar_cont_clst, tol = tol), center = T, scale = T),
    chloroquine_cont_scale_clst = scale(logit(chloroquine_cont_clst, tol = tol), center = T, scale = T),
    amodiaquine_cont_scale_clst = scale(logit(amodiaquine_cont_clst, tol = tol), center = T, scale = T),
    quinine_cont_scale_clst = scale(logit(quinine_cont_clst, tol = tol), center = T, scale = T),
    act_cont_scale_clst = scale(logit(act_cont_clst, tol = tol), center = T, scale = T),
    otherartm_cont_scale_clst = scale(logit(otherartm_cont_clst, tol = tol), center = T, scale = T),
    other_cont_scale_clst = scale(logit(other_cont_clst, tol = tol), center = T, scale = T),
    anyatm_cont_scale_clst = scale(logit(anyatm_cont_clst, tol = tol), center = T, scale = T)
    ) %>% 
  dplyr::rename(hv001 = v001)  %>% 
  dplyr::mutate(hv001 = as.numeric(hv001)) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))

dt <- dplyr::left_join(dt, kdsrv_fvr_clst, by = "hv001")

dt %>% 
  group_by(hv001) %>% 
  dplyr::summarise(n = sum(anyatm_cont_clst))

# some clusters missing kids with fevers, so have NAs

#.............
# Cluster Level Kid RDT/Micro
#.............
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ Percentage of Children Classified in Two Test
pr <- readRDS("data/raw_data/dhsdata/datasets/CDPR61FL.rds")
rdtmicro <- pr %>% 
  dplyr::filter(hc1 < 60 & hc1 >= 6) %>% # age
  dplyr::filter(hv042 == 1) %>% # household selected for Hb
  dplyr::filter(hv103 == 1) %>% #de facto
  dplyr::mutate(hv005_wi = hv005/1e6,
                hml32_fctb = haven::as_factor(hml32),
                hml32_fctb = if_else(hml32_fctb %in% c("positive", "negative"), 
                                     hml32_fctb, factor(NA)),
                hml32_fctb =  forcats::fct_drop(hml32_fctb),
                hml32_numb = ifelse(hml32_fctb == "positive", 1, ifelse(hml32_fctb == "negative", 0, NA)),
                
                hml35_fctb = haven::as_factor(hml35),
                hml35_fctb = if_else(hml35_fctb %in% c("positive", "negative"), 
                                     hml35_fctb, factor(NA)),
                hml35_fctb =  forcats::fct_drop(hml35_fctb),
                hml35_numb = ifelse(hml35_fctb == "positive", 1, ifelse(hml35_fctb == "negative", 0, NA))
  )


rdtmicro_srvy <- rdtmicro %>% srvyr::as_survey_design(ids = hv001, 
                                          strata = hv023, 
                                          weights = hv005_wi)

rdtmicro_sum <- rdtmicro_srvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    RDTprev_cont_clst = srvyr::survey_mean(hml35_numb, na.rm = T, vartype = c("se")),
    microprev_cont_clst = srvyr::survey_mean(hml32_numb, na.rm = T, vartype = c("se"))
  ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se"))) %>% 
  dplyr::mutate(RDTprev_cont_scale_clst = scale(logit(RDTprev_cont_clst, tol = tol)),
                microprev_cont_scale_clst = scale(logit(microprev_cont_clst, tol=tol)),
                hv001 = as.numeric(hv001))

dt <- dplyr::left_join(dt, rdtmicro_sum, by = "hv001")




#..........................................................................................
#                               Final Write Out
#..........................................................................................
keep <- colnames(dt)[ grepl("_cont|_fctm|_fctb|_scale", colnames(dt)) ]
keep <- c(keep, "hv001", "hv023", "hv005", "hv005_wi",
          "houseid", "hvdate_dtdmy", "hvyrmnth_dtmnth", "hvyrmnth_dtmnth_lag",
          "urban_rura", "latnum", "longnum", "geometry", "adm1name")


dtkp <- dt[, keep]

saveRDS(dtkp, file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")


