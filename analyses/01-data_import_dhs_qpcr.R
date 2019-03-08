#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data and merge PCR and DHS recodes for CD2013
# This script will focus on merging the DHS data and our PCR results
# Will then merge this to the geospatial and climate data in the 01-data_import_space file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
# DHS munging
devtools::install_github("OJWatson/rdhs", ref="master")
library(rdhs)


#---------------------------------------------------------------------------------
# Using rDHS to pull down CD2013
#---------------------------------------------------------------------------------
# https://ojwatson.github.io/rdhs/
rdhs::set_rdhs_config(email = "nbrazeau@med.unc.edu",
                      project = "Malaria Spatiotemporal Analysis",
                      config_path = "rdhs.json",
                      global = FALSE,
                      cache_path = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/")

survs <- rdhs::dhs_surveys(countryIds = c("CD"),
                           surveyYearStart = 2013)


datasets <- rdhs::dhs_datasets(surveyIds = survs$SurveyId,
                               fileFormat = "flat")

# download all DHS datsets 
downloads <- rdhs::get_datasets(datasets$FileName[!grepl("GC", datasets$FileName)]) # still not reading GC correctly



#---------------------------------------------------------------------------------
# Read in and merge qPCR data
#---------------------------------------------------------------------------------
# read in data
pfpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pf_alladults_v4.csv", 
                         col_names = T) %>% 
         magrittr::set_colnames(tolower(colnames(.))) %>% 
         dplyr::select(c("hivrecode_barcode", "pfldh", "fcq_mean")) %>% 
         dplyr::rename(pfctmean = fcq_mean)
  
  
pvpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pv_alladults_v2.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::select(c("hivrecode_barcode", "pv18s", "corrected_ct", "original_platemnum")) %>% 
  dplyr::rename(pvctcrrct = corrected_ct)

popcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Po_alladults_V2.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::select(c("hivrecode_barcode", "po18s", "corrected_ct_adj")) %>% 
  dplyr::rename(poctcrrct = corrected_ct_adj)

# bind together the pf, po, pv results
panplasmpcrres <- dplyr::inner_join(pfpcr, popcr, by="hivrecode_barcode") %>% 
  dplyr::inner_join(., pvpcr, by="hivrecode_barcode")

if(nrow(panplasmpcrres) != nrow(pfpcr) & nrow(pfpcr) != nrow(popcr) & nrow(popcr) != nrow(pvpcr)){
  stop("Abandon all faith, ye who enter here. Look upstream to how the CD2013 pcr calls were made. There has been an error in barcode renaming and/or merging")
}

#---------------------------------------------------------------------------------
# Read in and match PR barcode to qpcr data
#---------------------------------------------------------------------------------
pr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDPR61FL.rds")
ar <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDAR61FL.rds") %>% 
  dplyr::rename(hivrecode_barcode = hiv01) %>% 
  dplyr::mutate(hivrecode_barcode = gsub(" ", "", hivrecode_barcode),
                hivrecode_barcode = tolower(hivrecode_barcode))

# match HIV/PCR barcodes with the PR recode
# The HIV recoded barcodes are under the variables HA62 (females) in the PR recode
str(pr$ha62)
attr(pr$ha62, "label")
attr(pr$ha62, "labels")
# check for other odd values
pr$ha62[pr$ha62 %in% c("99991", as.character(seq(99997:99999)))]
head(xtabs(~pr$ha62))

# no barcode for the men
sum(factor(pr$ha62) == "     ")
sum(factor(pr$ha62) == "?    ")

# The HIV recoded barcodes are under the variables HB62 (males) in the PR recode
str(pr$hb62)
attr(pr$hb62, "label")
attr(pr$hb62, "labels")
# check for other odd values
pr$hb62[pr$hb62 %in% c("99991", as.character(seq(99997:99999)))]
head(xtabs(~pr$hb62))

# no barcode for the men
sum(factor(pr$hb62) == "     ")
sum(factor(pr$hb62) == "?    ")


# make merge table
infodf <- tibble(ha62 = c("     ",           "99991", "99992",     "99993",  "99994",       "99995",   "99996",  "?    ",   "     ", "     ",     "     ",  "     ",       "     ",   "     ",   "     ", "?    "),
                 hb62 = c("     ",           "     ", "     ",     "     ",  "     ",       "     ",   "     ",  "     ",   "99991", "99992",     "99993",  "99994",       "99995",   "99996",   "?    ", "?    "),
                 info = c("No info for F/M", "F unk", "F Incmplt", "F dmgd", "F not prsnt", "F rfsd",  "F othr", "F mssng", "M unk", "M Incmplt", "M dmgd", "M not prsnt", "M rfsd",  "M othr",  "M mssng", "both missing"))

# let's make the barcode column now and join with ar
arpr <- pr  %>% 
  dplyr::left_join(x=., y = infodf, by=c("ha62", "hb62")) %>%
  dplyr::mutate(hivrecode_barcode = ifelse(is.na(info) & hb62 != "     ", hb62,
                                 ifelse(is.na(info) & ha62 != "     ", ha62, NA)
                                 ),
                hivrecode_barcode = tolower(factor(hivrecode_barcode))
                ) %>% 
  inner_join(., ar, by = "hivrecode_barcode") 

# check to see if barcodes had missing in AR recode
if(TRUE %in% as.data.frame(table(arpr$info, arpr$hivrecode_barcode))[,3] != 0){
  stop("barcode parsing error with info missing and barcode")
}

# check to make sure every PR recode has an AR recode 
if(sum(!is.na(arpr$hivrecode_barcode)) != nrow(ar)){
  stop("predicted barcode number does not match barcode number in AR (HIV) recode")
}



#---------------------------------------------------------------------------------
# Read in MR and IR 
#---------------------------------------------------------------------------------
mr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDMR61FL.rds")
wr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDIR61FL.rds")

# make some new household variables for wrangling that are not in the PR but are important 
# note https://dhsprogram.com/data/Guide-to-DHS-Statistics/index.htm#t=HIV_Prevalence.htm&rhsearch=hiv&rhhlterm=hiv&rhsyns=%20
# The HIV test results data (AR file) should first be merged to the women’s dataset (IR file) or the men’s dataset (MR file) by cluster number (v001/mv001/hivclust), household number (v002/mv002/hivnumb) and line number (v003/mv003/hivline)

mrar <- ar %>% 
  dplyr::rename(mv001 = hivclust,
                mv002 = hivnumb,
                mv003 = hivline) %>% 
  dplyr::select(c("mv001", "mv002", "mv003", "hivrecode_barcode")) %>% 
  inner_join(., mr, by = c("mv001", "mv002", "mv003") )

wrar <- ar %>% 
  dplyr::rename(v001 = hivclust,
                v002 = hivnumb,
                v003 = hivline) %>% 
  dplyr::select(c("v001", "v002", "v003", "hivrecode_barcode")) %>% 
  inner_join(., wr, by = c("v001", "v002", "v003") )




 # WHAT IS GOING ON HERE?S??@??!
sum(is.na(mr[, c("mv001", "mv002", "mv003")]))
sum(is.na(wr[, c("v001", "v002", "v003")]))


nrow(mrar) + nrow(wrar) == nrow(ar)

potentialmatch <- as.data.frame(rbind(as.matrix(wr[, c("v001", "v002", "v003")]), 
                        as.matrix(mr[, c("mv001", "mv002", "mv003")])))
colnames(potentialmatch) <- c("hivclust", "hivnumb",  "hivline")

nomatch <- ar[ !ar$hiv01 %in% c(mrar$hiv01, wrar$hiv01), ]


mrprar <- left_join(x=arpr, y=mrar, by = "hivrecode_barcode")
wrmrprar <- left_join(x=mrprar, y=wrar, by = "hivrecode_barcode")


#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

gewrmrprar <- left_join(x=wrmrprar, y=ge, by = "hv001") 

#---------------------------------------------------------------------------------
# DHS Geospatial Covariates
#---------------------------------------------------------------------------------
gc <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/CDGC62FL/CDGC62FL.csv") %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

gcgewrmrprar <- dplyr::left_join(gewrmrprar, gc, by = c("dhsid", "dhscc", "dhsyear", "hv001"))



#---------------------------------------------------------------------------------
# FINAL -- subset the DHS data to just the PCR data 
#---------------------------------------------------------------------------------
panplasmpcrres_gcgewrmrprar <- dplyr::inner_join(panplasmpcrres, gcgewrmrprar, by = "hivrecode_barcode")




#---------------------------------------------------------------------------------
# write out
#---------------------------------------------------------------------------------

# write out joined HIV recode to PR, IR, MR, GE, and GC with our plasmodium PCR results 

saveRDS(panplasmpcrres_gcgewrmrprar, file = "~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")


# set up figure dir for later storage
if(!dir.exists("figures")){
  dir.create("figures")
}

