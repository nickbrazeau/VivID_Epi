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
                      cache_path = "/Users/NFB/Documents/GitHub/VivID_Epi")

survs <- rdhs::dhs_surveys(countryIds = c("CD"),
                           surveyYearStart = 2013)


datasets <- rdhs::dhs_datasets(surveyIds = survs$SurveyId,
                               fileFormat = "flat")

# download all DHS datsets 
downloads <- rdhs::get_datasets(datasets$FileName) 

#---------------------------------------------------------------------------------
# pull down maps
#---------------------------------------------------------------------------------
#spatial from GADM -- these are polygon files
drclvl0 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_0_sf.rds", httr::write_disk(path="data/gadm_drclvl0.rds", overwrite = T))
drclvl1 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_1_sf.rds", httr::write_disk(path="data/gadm_drclvl1.rds", overwrite = T))
drclvl2 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_2_sf.rds", httr::write_disk(path="data/gadm_drclvl2.rds", overwrite = T))
drclvl3 <- httr::GET(url = "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf/gadm36_COD_3_sf.rds", httr::write_disk(path="data/gadm_drclvl3.rds", overwrite = T))

DRCprov <- readRDS("data/gadm_drclvl1.rds")
colnames(DRCprov) <- tolower(colnames(DRCprov))
colnames(DRCprov)[4] <- "adm1name" # to match the DHS province names
# need to strip accent marks also to match the DHS province names
# https://stackoverflow.com/questions/20495598/replace-accented-characters-in-r-with-non-accented-counterpart-utf-8-encoding
# thanks to @Thomas for this great trick

unwanted_array = list(   'Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                         'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                         'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'ß'='Ss', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                         'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                         'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )

DRCprov$adm1name <- chartr(paste(names(unwanted_array), collapse=''),
       paste(unwanted_array, collapse=''),
       DRCprov$adm1name)

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


#---------------------------------------------------------------------------------
# Read in and merge qPCR data
#---------------------------------------------------------------------------------
# read in data
pfpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pf_alladults_v4.csv", 
                         col_names = T) %>% 
         magrittr::set_colnames(tolower(colnames(.))) %>% 
         dplyr::select(c("hivrecode_barcode", "pfldh"))
  
pvpcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Pv_alladults_v1.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::select(c("hivrecode_barcode", "pv18s"))

popcr <- readr::read_csv(file="/Volumes/share/1. Data/2. Data Set Processing/CD2013DHS_Adults_Construction/Po_alladults_V2.csv", 
                         col_names = T) %>% 
  magrittr::set_colnames(tolower(colnames(.))) %>% 
  dplyr::select(c("hivrecode_barcode", "po18s"))



# bind together the pf, po, pv results
panplasmpcrres <- dplyr::inner_join(pfpcr, popcr, by="hivrecode_barcode") %>% 
  dplyr::inner_join(., pvpcr, by="hivrecode_barcode")

if(nrow(panplasmpcrres) != nrow(pfpcr) & nrow(pfpcr) != nrow(popcr) & nrow(popcr) != nrow(pvpcr)){
  stop("Abandon all faith, ye who enter here. Look upstream to how the CD2013 pcr calls were made. There has been an error in barcode renaming and/or merging")
}

#---------------------------------------------------------------------------------
# Read in and match PR barcode to qpcr data
#---------------------------------------------------------------------------------
pr <- readRDS(file = "~/Documents/GitHub/VivID_Epi/datasets/CDPR61FL.rds")
ar <- readRDS(file = "~/Documents/GitHub/VivID_Epi/datasets/CDAR61FL.rds") %>% 
  dplyr::rename(hivrecode_barcode = hiv01) %>% 
  dplyr::mutate(hivrecode_barcode = gsub(" ", "", hivrecode_barcode))

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
                hivrecode_barcode = factor(hivrecode_barcode)
                ) %>% 
  inner_join(., ar, by = "hivrecode_barcode") %>%
  dplyr::mutate(hivrecode_barcode = tolower(hivrecode_barcode)
                )

# check to see if barcodes had missing in AR recode
if(TRUE %in% as.data.frame(table(arpr$info, arpr$hivrecode_barcode))[,3] != 0){
  stop("barcode parsing error with info missing and barcode")
}

# check to make sure every PR recode has an AR recode 
if(sum(!is.na(arpr$hivrecode_barcode)) != nrow(ar)){
  stop("predicted barcode number does not match barcode number in AR (HIV) recode")
}



# write out joined HIV recode to PR, can use this for panplasmodium results
if(!dir.exists(paths = "data")){
  dir.create("data")
}




save(DRCprov, ge, arpr, panplasmpcrres, file = "data/vividepi_raw.rda")


