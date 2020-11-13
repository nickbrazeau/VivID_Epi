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

# create bounding box of Central Africa for Speed/sanity checks
# https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952
caf <- as(raster::extent(10, 40,-18, 8), "SpatialPolygons")
sp::proj4string(caf) <- "+init=epsg:4326"

#--------------------------------------------------------------
# Section 1:Pulling data-map file for all recodes
#-------------------------------------------------------------- 
# CD2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# https://dhsprogram.com/data/Guide-to-DHS-Statistics/ -- note this is version 7
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")

# subset
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>%  # drop observations with missing geospatial data 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::filter(hv102 == 1) %>% # subset to de-jure https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm
  dplyr::filter(hiv05 != 0) # drop observations with samplings weights set to 0

# sanity check
sf::st_crs(dt)
# liftover to conform with rgdal updates http://rgdal.r-forge.r-project.org/articles/PROJ6_GDAL3.html
dt <- sp::spTransform(sf::as_Spatial(dt), CRSobj = sp::CRS("+init=epsg:4326"))
sp::identicalCRS(dt, caf)
# back to tidy 
dt <- sf::st_as_sf(dt)


#--------------------------------------------------------------
# Section 2: Looking at recodes, manual data wrangle
#-------------------------------------------------------------- 
#..........................
# Exposure of Interests
#..........................
# SURVEY CHARACTERISTICS/WEIGHTS
# 1. hiv05, hiv level weights
# 2. hv002, households


# COINFECTIONS/BIOMARKER VARIABLES
# 1. pfldh coinfection ; (personal)
# 2. HIV Status (HIV03) ; (personal)
# 3. Duffy phenotype (wetlab result -- nearly all Duffyneg, so won't be formally considered in model)

# SOCIOECOLOGICAL VARIABLES
# 1. Biological Sex (HV104)
# 2. Age (HV105)
# 3. Main floor material (categorical: HV213)
# 3. Main wall material (categorical: HV214)
# 3. Main roof material (categorical: HV215)
# 3 -> 4. Building material (recode)
# 5. Wealth Index (my recode; base wealth is hv270) 
# 6. Highest year of education completed (continous: HV108) 
# 7. Occupation (categorical: "hv717")
# 8. Number of Household Members (continuous: HV009)


# MALARIA-INTERVENTIONS
# 1. ITN Use
# 
# PHYSICAL/LANDSCAPE/CLIMATE VARIABLES
# 1. Cluster altitude (HV040)
# 2. Temparature (my recode)
# 3. Precipation (my recode)



#########################################################################################################
##############                             SURVEY CHARACTERISTICS                          ##############
#########################################################################################################
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hiv05_wi = hiv05/1e6)


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
  dplyr::mutate(hvdate_dtdmy = lubridate::dmy(paste(hv016, hv006, hv007, sep = "/")),
                hvyrmnth_dtmnth = paste(lubridate::year(hvdate_dtdmy), lubridate::month(hvdate_dtdmy), sep = "-"))


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
floor <- readr::read_csv("internal_datamap_files/pr_floor_liftover.csv")
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
tab <- xtabs(~dt$wlthrcde_fctm + haven::as_factor(dt$hv270)) # looks just OK. 
tab 
sum(diag(tab))/sum(tab) # 56% concordance. Presumambly, household type was a big driver given the 0 cells in the lower corner
# look at binary
xtabs(~dt$wlthrcde_fctb + haven::as_factor(dt$hv270))
xtabs(~dt$wlthrcde_fctb + dt$wlthrcde_fctm)

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
edu <- readr::read_csv("internal_datamap_files/pr_education_liftover.csv")
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
table(factor(haven::as_factor(dt$farmer_fctb)), useNA = "always")

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

#.............
# Weather
#.............
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
hlthdist_out <- readRDS("data/derived_data/hlthdist_out_wlk_trvltime.rds") 
# drop geom
sf::st_geometry(hlthdist_out) <- NULL 
# drop extra columns
hlthdist_out <- hlthdist_out %>% 
  dplyr::select(c("hv001", "hlthdist_cont_clst"))

dt <- dt %>% 
  dplyr::left_join(x=., y = hlthdist_out, by = "hv001") %>% 
  dplyr::mutate(
    hlthst_duration_cont_scale_clst = my.scale(hlthdist_cont_clst, center = T, scale = T),
    hlthst_duration_fctb_clst = ifelse(hlthdist_cont_clst > 60, "far", "near"),
    hlthst_duration_fctb_clst = factor(hlthst_duration_fctb_clst, levels = c("near", "far")))

# look at output
xtabs(~dt$hlthst_duration_fctb_clst + haven::as_factor(dt$hv270))
xtabs(~dt$hlthst_duration_fctb_clst + haven::as_factor(dt$urban_rura_fctb))

#..........................................................................................
#                               Final Write Out
#..........................................................................................
#...........................
# All Observations
#...........................
saveRDS(dt, file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")

#...........................
# Complete Observations
#...........................
# get final covariates
dcdr <- readxl::read_excel(path = "~/Documents/GitHub/VivID_Epi/model_datamaps/sub_DECODER_covariate_map_v3.xlsx", sheet = 1) %>% 
  dplyr::pull(c("column_name"))

dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
dt.cc <- dt  %>% 
  dplyr::select(c("pv18s", "pfldh", "po18s", dcdr)) %>% 
  dplyr::filter(complete.cases(.)) 

saveRDS(object = dt.cc,
        file = "~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode_completecases.rds")


