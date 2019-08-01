#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and clean the various covariates in
# the CD2013, weather, and other datasets
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
# 4. Anemia Level (HA57 & HB57)
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
# 1. ITN Use
# 2. Cluster level antimalarial use
# 
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)
# 3. Temparature (ecological import)
# 4. Precipation (ecological import)
# 6. Type of Residence i.e. Rural or Urban (recode)


#########################################################################################################
##############                             SURVEY CHARACTERISTICS                          ##############
#########################################################################################################
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6
  )

# #.............
# # dates
# #.............
# dt <- dt %>% 
#   dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")))
# 
# # NOTE, some clusters have survey start and end dates that are in two months 
# # (eg boundaries aren't clean/coinciding with a month. Naturally). Given
# # grouping by month, need to assign a clusters "month" on the majority of days 
# # that were spent surveying that clusters
# 
# # clusters without clean boundaries
# clst_mnth_bounds <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
#   dplyr::mutate(mnth = lubridate::month(hvdate_dtdmy)) %>% 
#   dplyr::group_by(hv001) %>% 
#   dplyr::summarise(moremnths = length(unique(mnth))) %>% 
#   dplyr::filter(moremnths > 1)
# 
# clst_mnth_bounds.assign <- dt[, c("hv001", "hvdate_dtdmy")] %>% 
#   dplyr::filter(hv001 %in% clst_mnth_bounds$hv001) %>% 
#   dplyr::mutate(hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-")) %>% 
#   dplyr::group_by(hv001, hvyrmnth_dtmnth) %>% 
#   dplyr::summarise(n = n()) %>% 
#   dplyr::filter(n == max(n)) %>% 
#   dplyr::select(-c("n"))
# 
# sf::st_geometry(clst_mnth_bounds.assign) <- NULL
# 
# dt <- dt %>% 
#   dplyr::left_join(x=., y = clst_mnth_bounds.assign, by = "hv001") %>% 
#   dplyr::mutate(hvyrmnth_dtmnth = ifelse(is.na(hvyrmnth_dtmnth),
#                                          paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"),
#                                          hvyrmnth_dtmnth))
# 
# dates <- readr::read_csv("internal_datamap_files/pr_date_liftover.csv")
# dt <- dt %>% 
#   dplyr::left_join(x=., y=dates, by = "hvyrmnth_dtmnth") %>% 
#   dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
# 
# xtabs(~dt$hvyrmnth_dtmnth + dt$hvyrmnth_dtmnth_lag)


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
    po18sct_cont = poctcrrct,
    po18sct_cont_log = log(poctcrrct + tol),
    pv18sct_cont = pvctcrrct,
    pv18sct_cont_log = log(pvctcrrct + tol),
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
# Anemia/Hemoglobin
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
                hab57_fctb = factor(hab57_fctb, levels = c("no", "yes"))
                )
            
# check recode
xtabs(~dt$hab57_fctb + haven::as_factor(dt$ha57) + haven::as_factor(dt$hv104), addNA = T)  
xtabs(~dt$hab57_fctb + haven::as_factor(dt$hb57) + haven::as_factor(dt$hv104), addNA = T)  

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

wallfloorroofrecode <- dt %>% 
  dplyr::select(c("hv213_liftover", "hv214_liftover", "hv215_liftover"))
sf::st_geometry(wallfloorroofrecode) <- NULL
dt <- dt %>% 
  dplyr::mutate(
    housecount = apply(wallfloorroofrecode, 1, function(x){return(sum(x == "non-rudimentary"))}),
    hv21345_fctb = ifelse(housecount == 3, "modern", "traditional"), # per PMID: 28222094
    hv21345_fctb = factor(hv21345_fctb))
                
# check -- seems reasonable
xtabs(~dt$hv21345_fctb, addNA = T)


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
xtabs(~dt$wlthrcde_fctb + haven::as_factor(dt$hv270)) # looks OK. Some poor/poorest got placed higher than expected



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
  dplyr::mutate( hv106_fctb = factor(hv106_fctb, levels = c("lower", "higher")) )

xtabs(~dt$hv106 + dt$hv106_fctb, addNA = T)

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
# ethnic <-  readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_ethnic_liftover.csv")
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

# #------------------------------------------
# # Occupation
# #------------------------------------------
# # Going to dichotomize as suspected outdoor versus indoor
# # Will code "others" as NA because cannot discern, only 2 obs
# 
# dt$hab717 <- ifelse(haven::as_factor(dt$hv104) == "female", dt$v717, dt$mv717)
# table(dt$v717)
# table(dt$mv717)
# table(dt$hab717)
# 
# occupation <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_occupation_liftover.csv")
# occupation$hab717 <- factor(occupation$hv717) # my mac being a pain
# 
# dt <- dt %>% 
#   dplyr::mutate(hab717 = haven::as_factor(hab717), 
#                 hab717 =  forcats::fct_drop(hab717))
# dt <- dt %>%
#   left_join(x=., y=occupation, by="hab717") %>% 
#   dplyr::mutate(hab717_fctb = factor(jobconditions),
#                 hab717_fctb = forcats::fct_relevel(hab717_fctb, "indoors"))
# 
# xtabs(~dt$hab717 + dt$hab717_fctb, addNA = T)

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
# #.............
# # Health Insurance
# #.............
# dt <- dt %>% 
#   dplyr::mutate(hab481 =  ifelse(haven::as_factor(hv104) == "female", v481, mv481),
#                 hab481 = ifelse(hab481 == 9, NA, hab481),
#                 hab481_fctb = factor(hab481, levels = c(0,1), labels = c("no", "yes")))
# xtabs(~dt$hab481_fctb, addNA = T)

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
      !(haven::as_factor(hml4) %in% c("36", "more than 3 years ago", "don't know", "missing")) & haven::as_factor(hml20) == "yes", 1,
      # (ii) conventional ITNs that were 1 y old 
      ifelse(haven::as_factor(hml4) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1, 
             # or were retreated within the year before the survey
             ifelse(haven::as_factor(hml9) %in% c(0:12) & haven::as_factor(hml12) == "only treated (itn) nets", 1, 
                    0))), # we know no missing net from above. they either reported yes or no at some level
    ITN_fctb = factor(ITN_fctb, levels = c(0,1), labels = c("no", "yes")),
    ITN_fctb = forcats::fct_drop(ITN_fctb)
  )


# sanity check
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml10)) 
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml12)) 
xtabs(~dt$ITN_fctb + haven::as_factor(dt$hml20)) 

# #.............
# # LLIN-type of Inseciticide for INDIVIDUAL
# # Note, must have LLIN to have insecticide (120 missing LLIN insecticide types, 8500 no LLIN)
# #.............
# # read insecticide liftover table
# insctcd <- readr::read_csv("~/Documents/GitHub/VivID_Epi/internal_datamap_files/pr_insecticide_liftover.csv")
# 
# dt <- dt %>%
#   dplyr::mutate(hml7 = haven::as_factor(hml7)) %>%
#   left_join(x=., y=insctcd, by="hml7") %>%
#   dplyr::mutate(insctcd_fctm = factor(ifelse(haven::as_factor(hml20) == "no", "none", insctcd)),
#                 insctcd_fctm = forcats::fct_relevel(insctcd_fctm, "none")
#   )
# 
# 
# # sanity checks
# xtabs(~dt$insctcd_fctm + haven::as_factor(dt$hml20), addNA=T)
# xtabs(~dt$insctcd_fctm + dt$hml7, addNA=T)
# xtabs(~dt$hml20_fctb + dt$hml7, addNA=T)
# 
# 
# # is this confounded by age?
# xtabs(~dt$hml4 + dt$insctcd_fctm, addNA=T) # nope



#########################################################################################################
##############                             CLUSTER LEVEL VARIABLES                         ##############
#########################################################################################################
# ALL CLUSTER LEVEL VARIABLES WILL BE APPENDED WTIH "_clst"  
dtsrvy <- makecd2013survey(survey = dt)

#..........................................................................................
#                                  COINFECTIONS/BIOMARKER VARIABLES
#..........................................................................................
# Hbmiss <- dt[is.na(dt$hab56_cont),]
# table(Hbmiss$hv001) # looks well dispersed among clusters
# 
# coinfxnbiomrk <- dtsrvy %>% 
#   dplyr::group_by(hv001) %>% # cluster level 
#   dplyr::summarise(
#     pfldh_cont_clst = srvyr::survey_mean(x = pfldh),
#     po18s_cont_clst = srvyr::survey_mean(x = po18s),
#     hiv03_cont_clst = srvyr::survey_mean(x = hiv03),
#     hab56_cont_clst = srvyr::survey_mean(hab56_cont, na.rm = T)) %>% 
#   dplyr::mutate(
#     pfldh_cont_scale_clst = my.scale(pfldh_cont_clst, center = T, scale = T),
#     po18s_cont_scale_clst = my.scale(po18s_cont_clst, center = T, scale = T),
#     hiv03_cont_scale_clst = my.scale(hiv03_cont_clst, center = T, scale = T),
#     hab56_cont_scale_clst = my.scale(hab56_cont_clst, center = T, scale = T)
#   ) %>% 
#   dplyr::select(-c(dplyr::ends_with("_se")))
# 
# sapply(coinfxnbiomrk, summary) # looks clean, had to remove NAs in Hb because 49 missing (spread across 39 clusters)
# 
# dt <- dplyr::left_join(x = dt, y = coinfxnbiomrk)


#..........................................................................................
#                                ECOLOGICAL VARIABLES
#..........................................................................................
source("R/00-functions_maps.R")
# #.............
# # Precipitation (CHRIPS) & Temperature (Emch/Manny)
# #.............
# # precip data
# precip <- list.files(path = "data/raw_data/weather_data/CHIRPS/", full.names = T)
# precipfrst <- lapply(precip, readRasterBB, bb = bb)
# precipdf <- tibble::tibble(names = basename(precip)) %>% 
#   dplyr::mutate(names = gsub("chirps-v2.0.|.tif", "", names),
#                 year = stringr::str_split_fixed(names, "\\.", n=2)[,1] ,
#                 mnth =  stringr::str_split_fixed(names, "\\.", n=2)[,2] ,
#                 hvdate_dtdmy = lubridate::dmy(paste(1, mnth, year, sep = "/")),
#                 year = lubridate::year(hvdate_dtdmy),
#                 mnth = lubridate::month(hvdate_dtdmy),
#                 hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
#                 precipraster = precipfrst) %>% 
#   dplyr::select(c("hvyrmnth_dtmnth_lag", "precipraster"))
# 
# # temp data
# temp <- raster::stack("data/raw_data/weather_data/emch_manny/cru_ts4.01.2011.2016.tmp.dat.nc") %>% 
#         raster::crop(x = ., y = sf::as_Spatial(bb))
# 
# tempdf <- tibble::tibble(orignnames = names(temp),
#                          names = gsub("X", "", names(temp)),
#                          hvdate_dtdmy = lubridate::ymd(names),
#                          year = lubridate::year(hvdate_dtdmy),
#                          mnth = lubridate::month(hvdate_dtdmy),
#                          hvyrmnth_dtmnth_lag = factor(paste(year, mnth, sep = "-")),
#                          tempraster = lapply(orignnames, raster::subset, x = temp)) %>% 
#   dplyr::select(c("hvyrmnth_dtmnth_lag", "tempraster"))
# 
# 
# 
# wthrnd <- dt[,c("hv001", "hvyrmnth_dtmnth_lag", "geometry", "urban_rura")] %>% 
#   dplyr::mutate(buffer = ifelse(urban_rura == "R", 10, 2))
# wthrnd <- wthrnd[!duplicated(wthrnd$hv001),]
# 
# wthrnd <- wthrnd %>% 
#   dplyr::left_join(., tempdf) %>% 
#   dplyr::left_join(., precipdf)
# 
# # Drop in a for loop
# wthrnd$precip_lag_cont_clst <- NA
# wthrnd$temp_lag_cont_clst <- NA
# 
# for(i in 1:nrow(wthrnd)){
#   # precip
#   wthrnd$precip_lag_cont_clst[i] <- 
#     raster::extract(x = wthrnd$precipraster[[i]],
#                     y = sf::as_Spatial(wthrnd$geometry[i]),
#                     buffer = wthrnd$buffer[i],
#                     fun = mean,
#                     sp = F
#                     )
#   
#   # temp
#   wthrnd$temp_lag_cont_clst[i] <- 
#     raster::extract(x = wthrnd$tempraster[[i]],
#                     y = sf::as_Spatial(wthrnd$geometry[i]),
#                     buffer = wthrnd$buffer[i],
#                     fun = mean,
#                     sp = F
#     )
#   
# }
# 
# wthrnd <- wthrnd %>% 
#   dplyr::select(c("hv001", "hvyrmnth_dtmnth_lag", "precip_lag_cont_clst", "temp_lag_cont_clst")) %>% 
#   dplyr::mutate(hvyrmnth_dtmnth_lag = factor(hvyrmnth_dtmnth_lag))
# sf::st_geometry(wthrnd) <- NULL
# dt <- dt %>% 
#   dplyr::left_join(., wthrnd, by = c("hv001", "hvyrmnth_dtmnth_lag")) %>% 
#   dplyr::mutate(
#                 precip_lag_cont_scale_clst = my.scale(precip_lag_cont_clst, center = T, scale = T),
#                 temp_lag_cont_scale_clst = my.scale(temp_lag_cont_clst, center = T, scale = T)
#                 )
# 
# 

dt <- dt %>% 
  dplyr::mutate(precip_ann_cont_clst = ifelse(annual_precipitation_2015 == 9999, NA, annual_precipitation_2015), # note no missing (likely dropped with missing gps)
                precip_ann_cont_scale_clst = my.scale(precip_ann_cont_clst, center = T, scale = T),
                temp_ann_cont_clst = ifelse(night_land_surface_temp_2015 == 9999, NA, night_land_surface_temp_2015), # note no missing (likely dropped with missing gps)
                temp_ann_cont_scale_clst = my.scale(night_land_surface_temp_2015, center = T, scale = T)
  )


#.............
# Cluster-Level Altitude
#.............
dt <- dt %>% 
  dplyr::mutate(alt_dem_cont_clst = ifelse(alt_dem == 9999, NA, alt_dem), # note no missing (likely dropped with missing gps)
                alt_dem_cont_scale_clst = my.scale(alt_dem_cont_clst, center = T, scale = T)
  )

# #.............
# # Ape Habitat Overlap
# #.............
# drc_ape <- readRDS(file = "data/redlist_species_data_primate/drc_ape.rds")
# drc_ape_simplified <- rgeos::gSimplify(
#   sf::as_Spatial(sf::st_union(drc_ape)), 
#   tol = 1e-3)
# 
# dt <- dt %>% 
#   dplyr::mutate(ape_habitat_fctb_clst = factor(
#     ifelse(rgeos::gContains(spgeom1 = drc_ape_simplified,
#                             spgeom2 = sf::as_Spatial(sf::st_as_sf(.)),
#                             byid = T),
#            1, 0),
#     levels = c(0,1), labels = c("NoOverlap", "Overlap"))
#   )
# xtabs(~dt$ape_habitat_fctb_clst)

#.............
# Distance to Water Source
#.............
wtrdist_out <- readRDS("data/derived_data/hotosm_waterways_dist.rds")
dt <- dt %>% 
  dplyr::left_join(x=., y = wtrdist_out, by = "hv001") %>% 
  dplyr::mutate(wtrdist_cont_scale_clst = my.scale(wtrdist_cont_clst, center = T, scale = T)
                )


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
xtabs(~dt$urbanscore_fctm_clst + haven::as_factor(dt$hv025), addNA = T) # looks OK. Some rural places are now urban, which is consistent with what I expected

#.............
# Distance to Health Site
#.............
hlthdist_out <- readRDS("data/derived_data/hotosm_healthsites_dist.rds")
dt <- dt %>% 
  dplyr::left_join(x=., y = hlthdist_out, by = "hv001") %>% 
  dplyr::mutate(hlthdist_cont_scale_clst = my.scale(log(hlthdist_cont_clst + tol), center = T, scale = T)
                )

#.............
# ITN Cluster Usage; mean wealth; prop educ
#.............
summary(dt$hml20) # no missing
summary(dt$wlthrcde_fctm)
summary(dt$hv106_fctb)

democlust <- dtsrvy %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    ITN_cont_clst = srvyr::survey_mean((ITN_fctb == "yes"), vartype = c("se")),
    hab57_cont_clst = srvyr::survey_mean((hab57_fctb == "yes"), vartype = c("se")),
    hv106_cont_clst = srvyr::survey_mean((hv106_fctb == "higher"), vartype = c("se"), na.rm = T)
  ) %>% 
  dplyr::mutate(
    hml20_cont_scale_clst = my.scale(logit(ITN_cont_clst, tol = tol), center = T, scale = T),
    hv106_cont_scale_clst = my.scale(logit(hv106_cont_clst, tol = tol), center = T, scale = T)
          ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))

# find median wealth and then recode it as a factor 
clstwlth <- dtsrvy %>% 
  dplyr::group_by(hv001) %>%
  dplyr::summarise(
    wlthrcde_fctm_clst = srvyr::survey_quantile(as.numeric(wlthrcde_fctm), quantiles = 0.5)
  ) %>% 
  dplyr::select(c("hv001", "wlthrcde_fctm_clst_q50")) %>% 
  dplyr::mutate(
    wlthrcde_fctm_clst_q50 = floor(wlthrcde_fctm_clst_q50) # in case of ties
  )

# decode this
wlthliftover <- data.frame(
  wlthrcde_fctm_clst_q50 = 1:5,
  wlthrcde_fctm_clst = levels(dt$wlthrcde_fctm)
)

clstwlth <- dplyr::left_join(clstwlth, wlthliftover, by = "wlthrcde_fctm_clst_q50") %>% 
  dplyr::select(-c("wlthrcde_fctm_clst_q50"))

# now merge back into democlust
democlust <- dplyr::left_join(democlust, clstwlth, by = "hv001")


sapply(democlust, summary) # looks clean

dt <- dplyr::left_join(x = dt, y = democlust)


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
    actuse.scale = mean(anyatm_cont_scale_clst)) # looks good



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
  dplyr::mutate(RDTprev_cont_scale_clst = my.scale(logit(RDTprev_cont_clst, tol = tol)),
                microprev_cont_scale_clst = my.scale(logit(microprev_cont_clst, tol=tol)),
                hv001 = as.numeric(hv001))

dt <- dplyr::left_join(dt, rdtmicro_sum, by = "hv001")




#..........................................................................................
#                               Final Write Out
#..........................................................................................
# keep <- colnames(dt)[ grepl("_cont|_fctm|_fctb|_scale", colnames(dt)) ]
# keep <- c("pv18s", keep, "hv001", "hv023", "hv005", "hv005_wi",
#           "houseid", "hvdate_dtdmy", "hvyrmnth_dtmnth", "hvyrmnth_dtmnth_lag",
#           "urban_rura", "latnum", "longnum", "geometry", "adm1name")
# 
# 
# dtkp <- dt[, keep]

saveRDS(dt, file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")


