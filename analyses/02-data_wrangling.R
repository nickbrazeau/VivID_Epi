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
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3
set.seed(48)


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

dt[, colnames(dt)[grepl("pv18s|pfldh", colnames(dt))] ] %>% 
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
# wet lab. basically all individuals except for 3 -- no formal modeling needed (since no variation)


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

#.............
# Temperature 
#.............
dt <- dt %>% 
  dplyr::mutate(temp_ann_cont_clst = ifelse(night_land_surface_temp_2015 == 9999, NA, night_land_surface_temp_2015), # note no missing (likely dropped with missing gps)
                temp_ann_cont_scale_clst = my.scale(temp_ann_cont_clst, center = T, scale = T))
#.............
# Precipitation 
#.............
prcp <- readRDS("data/derived_data/annual_precipitation_2015_imputed.RDS")
dt <- dt %>% 
  dplyr::left_join(., prcp, by = "hv001") %>% 
  dplyr::rename(precip_ann_cont_clst = annual_precipitation_2015imp) %>% 
  dplyr::mutate(
    precip_ann_cont_scale_clst = my.scale(precip_ann_cont_clst, center = T, scale = T)
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
# going to use those four variables as latent variables of urbanicity


# urb <- readRDS(file = "data/derived_data/vividepi_urban_recoded.rds") %>% 
#   dplyr::select(c("hv001", "urbanscore")) %>% 
#   dplyr::rename(urbanscore_cont_clst = urbanscore)
# dt <- dplyr::left_join(dt, urb, by = "hv001")
# boxplot(dt$urbanscore ~ haven::as_factor(dt$hv025)) # looks OK. Some urban places look pretty rural... which is more or less what I expected

urb <- readRDS(file = "data/derived_data/vividepi_urban_recoded.rds")
dt <- dt %>% 
  dplyr::left_join(x=., y = urb, by = c("hv001", "hv025")) 


#.............
# Distance to Health Site
#.............
hlthdist_out <- readRDS("data/derived_data/hotosm_healthsites_dist.rds")
dt <- dt %>% 
  dplyr::left_join(x=., y = hlthdist_out, by = "hv001") %>% 
  dplyr::mutate(hlthdist_cont_scale_clst = my.scale(log(hlthdist_cont_clst + tol), center = T, scale = T)
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
    actuse.scale = mean(anyatm_cont_scale_clst)) # looks good




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


