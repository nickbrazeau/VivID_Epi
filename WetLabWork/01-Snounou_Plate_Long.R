#---------------------------------------------------------------------------
# Purpose of this script is to take Snounou plates to long format
# will also randomly select 10% for HRM confirmation by Sanger Sequencing
#---------------------------------------------------------------------------

#---------------
# Dependencies and imports
#---------------
library(tidyverse)
library(zoo)
set.seed(44)
sn <- readxl::read_excel(path =  "/Volumes/share/1. Data/1. Raw Data/Adult_Pv18s/VivID_Snounou_plate_maps.xlsx", sheet = 2, col_names = F) # Cedar/Nick's results Pv Snounou Plate Maps 

#---------------
# Plate Map Spread
#---------------
sn <- sn %>% 
  magrittr::set_colnames(c("plate", "rowlett", seq(1:12), "drop", "drop2")) %>% 
  dplyr::select(-c(drop, drop2)) %>% 
  dplyr::mutate(plate = zoo::na.locf(plate)) %>% 
  dplyr::filter(!is.na(rowlett)) %>% # drop column numbering for 96 well plates, i.e. where ther is no row letter
  tidyr::gather(key = "columnnum", "longbarcode", 3:ncol(.)) %>% # reshape
  dplyr::filter(!is.na(longbarcode)) %>%  # drop rows that don't have barcodes
  dplyr::mutate(original_plate = stringr::str_split_fixed(longbarcode, "_", n=3)[,1],
                original_plate_loc = stringr::str_split_fixed(longbarcode, "_", n=3)[,2],
                barcode = stringr::str_split_fixed(longbarcode, "_", n=3)[,3]) %>% 
  dplyr::mutate(barcode = tolower(barcode)) # all lower cases, to remove Case Sensitivity



adult_snounou <- readxl::read_excel(path = "/Volumes/share/1. Data/1. Raw Data/Adult_Pv18s/Snounou_Verification_v3.xlsx", sheet = 1)
adult_snounou <- adult_snounou %>% 
  dplyr::mutate(original_plate = stringr::str_split_fixed(SampleName, "_", n=3)[,1],
                              original_plate_loc = stringr::str_split_fixed(SampleName, "_", n=3)[,2],
                              barcode = stringr::str_split_fixed(SampleName, "_", n=3)[,3]) %>% 
  dplyr::mutate(barcode = tolower(barcode)) # all lower cases, to remove Case Sensitivity

sn$barcode[! sn$barcode %in% adult_snounou$barcode ]

#---------------
# Randomly select 10% of samples actually being used in this study from the 16,363
#---------------
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
pvpos <- dt[dt$pv18s_fctb == "viv+",]
pvbarcodes <- sample(pvpos$hivrecode_barcode, size = 0.1*nrow(pvpos), replace = F) %>% 
  tibble(barcode = .)


#---------------
# Find matches in SN
#---------------
sangerduffy <- inner_join(x=sn, y=pvbarcodes, by = "barcode")
if(nrow(sangerduffy) != floor(0.1*nrow(pvpos))){
  stop("There was a barcode matching error. See which one it has a one-off base with as you did in the original data generation")
}

sangerduffy <- sangerduffy %>% 
  dplyr::arrange(plate, rowlett, columnnum) %>% 
  dplyr::select(c("plate", "rowlett", "columnnum", "barcode")) %>% 
  write.csv(., file = "WetLabWork/confirmatory_duffy_sanger_smpls.csv",
            quote = F, row.names = F)


