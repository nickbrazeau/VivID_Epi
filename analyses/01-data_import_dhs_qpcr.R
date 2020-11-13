#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data and merge PCR and DHS recodes for CD2013
# This script will focus on merging the DHS data and our PCR results
# Will then merge this to the geospatial and climate data in the 01-data_import_space file
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
# DHS wrangling
remotes::install_github("OJWatson/rdhs", ref="master")
library(rdhs)


#---------------------------------------------------------------------------------
# Using rDHS to pull down CD2013
#---------------------------------------------------------------------------------
# https://ojwatson.github.io/rdhs/
rdhs::set_rdhs_config(email = "nbrazeau@med.unc.edu",
                      project = "Malaria Spatiotemporal Analysis",
                      config_path = "rdhs.json",
                      global = FALSE,
                      cache_path = "data/raw_data/dhsdata/")

survs <- rdhs::dhs_surveys(countryIds = c("CD"),
                           surveyYearStart = 2013)


datasets <- rdhs::dhs_datasets(surveyIds = survs$SurveyId,
                               fileFormat = "flat")

# download all DHS datsets 
downloads <- rdhs::get_datasets(datasets$FileName) 



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


####################################################################################
##########                          DHS MERGING                           ##########  
####################################################################################
# will follow this
#  https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm

#---------------------------------------------------------------------------------
# Read in and match PR barcode to qpcr data
#---------------------------------------------------------------------------------
pr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDPR61FL.rds")
ar <- readRDS(file = "data/raw_data/dhsdata/datasets/CDAR61FL.rds") 
class(ar) <- "data.frame" # drop this custom dhs class so it doesn't interfere w/ tidy
ar <- ar %>% 
  dplyr::rename(hv001 = hivclust,
                hv002 = hivnumb,
                hvidx = hivline,
                hivrecode_barcode = hiv01) %>% 
  dplyr::mutate(hivrecode_barcode = gsub(" ", "", hivrecode_barcode),
                hivrecode_barcode = tolower(hivrecode_barcode)) # rename and fix barcode for ar
# match HIV/PCR barcodes with the PR recode
arpr <- dplyr::inner_join(ar,pr, by = c("hv001", "hv002", "hvidx"))
arpr <- arpr %>% 
  dplyr::mutate(hmid = paste(hv001, hv002, hvidx, sep = "_"))

#---------------------------------------------------------------------------------
# Read in MR and IR 
#---------------------------------------------------------------------------------
mr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDMR61FL.rds")
class(mr) <- "data.frame"
mr <- mr %>% 
  dplyr::rename(hv001 = mv001,
                hv002 = mv002,
                hvidx = mv003,
                caseid = mcaseid) %>% 
  dplyr::select(c("hv001", "hv002", "hvidx", "caseid", "mv717"))

wr <- readRDS(file = "data/raw_data/dhsdata/datasets/CDIR61FL.rds") 
class(wr) <- "data.frame"
wr <- wr %>% 
  dplyr::rename(hv001 = v001,
                hv002 = v002,
                hvidx = v003)  %>% 
  dplyr::select(c("hv001", "hv002", "hvidx", "caseid", "v717"))

mrwr <- dplyr::full_join(mr, wr)
mrwr <- mrwr %>%
  dplyr::mutate(hmid = paste(hv001, hv002, hvidx, sep = "_"))

wrmrprar <- dplyr::left_join(arpr, mrwr)


#---------------------------------------------------------------------------------
# DHS spatial shapes/data
#---------------------------------------------------------------------------------

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::rename(hv001 = dhsclust) %>%  # for easier merge with PR
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) 

saveRDS(object = ge, file = "data/raw_data/dhsdata/VivIDge.RDS")


gewrmrprar <- dplyr::left_join(x=wrmrprar, y=ge, by = "hv001") 

#---------------------------------------------------------------------------------
# DHS Geospatial Covariates
#---------------------------------------------------------------------------------
gc <- readRDS("data/raw_data/dhsdata/datasets/CDGC62FL.rds") %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR

gcgegewrmrprar <- dplyr::left_join(gewrmrprar, gc, by = c("dhsid", "dhscc", "dhsyear", "hv001"))



#---------------------------------------------------------------------------------
# FINAL -- subset the DHS data to just the PCR data 
#---------------------------------------------------------------------------------
panplasmpcrres_gcgeprar <- dplyr::inner_join(panplasmpcrres, gcgegewrmrprar, by = "hivrecode_barcode")




#---------------------------------------------------------------------------------
# write out
#---------------------------------------------------------------------------------

# write out joined HIV recode to PR, IR, MR, GE, and GC with our plasmodium PCR results 
# in order to perserve sf features, covert to sf (and dataframe)
panplasmpcrres_gcgeprar <- sf::st_as_sf(panplasmpcrres_gcgeprar)
saveRDS(panplasmpcrres_gcgeprar, file = "data/raw_data/vividpcr_dhs_raw.rds")




