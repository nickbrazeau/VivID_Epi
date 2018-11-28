#----------------------------------------------------------------------------------------------------
# Purpose of this script is to import data and merge PCR and DHS recodes for CD2013
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
                      cache_path = "/Users/nickbrazeau/Documents/GitHub/VivID_Epi")

survs <- dhs_surveys(countryIds = c("CD"),
                     surveyYearStart = 2013)


datasets <- dhs_datasets(surveyIds = survs$SurveyId,
                         fileFormat = "flat")
# download all DHS datsets 
downloads <- get_datasets(datasets$FileName) 

#---------------------------------------------------------------------------------
# Read in qPCR data
#---------------------------------------------------------------------------------
# read in data
pfpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pf_alladults_v4.csv", 
                         col_names = T) %>% 
         magrittr::set_colnames(tolower(colnames(.))) %>% 
         dplyr::rename(barcode = hivrecode_barcode) %>% 
         dplyr::select(c("barcode", "pfldh"))
  
pvpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pv_alladults_v1.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(barcode = hivrecode_barcode) %>% 
  dplyr::select(c("barcode", "pv18s"))

popcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Po_alladults_V2.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::rename(barcode = hivrecode_barcode) %>% 
  dplyr::select(c("barcode", "po18s"))

















cd2013dhspr <- readRDS(file = "~/Documents/GitHub/CD2013DHS_Adults_qPCR_Curation/Pfalciparum/datasets/CDPR61FL.rds")
cd2013dhsar <- readRDS(file = "~/Documents/GitHub/CD2013DHS_Adults_qPCR_Curation/Pfalciparum/datasets/CDAR61FL.rds")
colnames(cd2013dhsar) <- c("cluster", "household", "Line", "barcode", "HIV02", "HIV03", "HIV05")



# match HIV/PCR barcodes with the PR recode
# The HIV recoded barcodes are under the variables HA62 (males) and HB62 (females)
str(cd2013dhspr$ha62)
attr(cd2013dhspr$ha62, "label")
attr(cd2013dhspr$ha62, "labels")
# of note, there are values in the ha62 and hb62 barcodes that are not listed. Assuming it is missing/damaged.
cd2013dhspr$ha62[cd2013dhspr$ha62 %in% c("99991", as.character(seq(99997:99999)))]
head(xtabs(~cd2013dhspr$ha62))
# The HIV recoded barcodes are under the variables HA62 (males) and HB62 (females)
str(cd2013dhspr$hb62)
attr(cd2013dhspr$hb62, "label")
attr(cd2013dhspr$hb62, "labels")
# of note, there are values in the ha62 and hb62 barcodes that are not listed. Assuming it is missing/damaged, etc.
cd2013dhspr$hb62[cd2013dhspr$hb62 %in% c("99991", as.character(seq(99997:99999)))]
head(xtabs(~cd2013dhspr$hb62))

# no barcode for the men
sum(factor(cd2013dhspr$hb62) == "     ")
sum(factor(cd2013dhspr$hb62) == "?    ")


# make merge table
infodf <- tibble(ha62 = c("     ",           "99991", "99992",     "99993",  "99994",       "99995",   "99996",  "?    ",   "     ", "     ",     "     ",  "     ",       "     ",   "     ",   "     "),
                 hb62 = c("     ",           "     ", "     ",     "     ",  "     ",       "     ",   "     ",  "     ",   "99991", "99992",     "99993",  "99994",       "99995",   "99996",   "?    "),
                 info = c("No info for F/M", "F unk", "F Incmplt", "F dmgd", "F not prsnt", "F rfsd",  "F othr", "F mssng", "M unk", "M Incmplt", "M dmgd", "M not prsnt", "M rfsd",  "M othr",  "M mssng"))

# let's make the barcode column now
cd2013dhspr <- cd2013dhspr %>%
  dplyr::left_join(x=cd2013dhspr, y = infodf, by=c("ha62", "hb62")) %>%
  dplyr::mutate(barcode = ifelse(is.na(info) & hb62 != "     ", hb62,
                                 ifelse(is.na(info) & ha62 != "     ", ha62, NA)
                                 )
                )

# check
if(TRUE %in% as.data.frame(table(cd2013dhspr$info, cd2013dhspr$barcode))[,3] != 0){
  stop("barcode parsing error with info missing and barcode")
}
if(sum(!is.na(cd2013dhspr$barcode)) != nrow(cd2013dhsar)){
  stop("predicted barcode number does not match barcode number in AR (HIV) recode")
}

# filter down to barcodes that are meaningful
colnames(cd2013dhspcr)[colnames(cd2013dhspcr) == "HIVrecode_barcode"] <- "barcode"
cd2013dhspr$barcode <- tolower(cd2013dhspr$barcode)
cd2013dhspr <- inner_join(cd2013dhspcr, cd2013dhspr, by=c("barcode"))

# add in GPS points
cd2013gps <- tibble(cluster = cd2013dhsge$DHSCLUST, lat = cd2013dhsge$LATNUM, long = cd2013dhsge$LONGNUM)
# Clusters with Lat and Long of 0,0 were not able to be identified and should have coordinates set to NA
cd2013gps <- cd2013gps %>%
  dplyr::mutate(long = ifelse(long == 0 & lat == 0, NA, long)) %>%
  dplyr::mutate(lat = ifelse(long == 0 & lat == 0, NA, lat))

cd2013dhspr <- left_join(x = cd2013dhspr, y=cd2013gps, by=c("cluster"))




save(cd2013dhspr, cd2013dhsge, drclong, file = "data/EPBID_raw.rda")




