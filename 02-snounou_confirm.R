#---------------------------------------------------------------------------
# Purpose of this script is to randomly select 94 samples
# (e.g. one plate worth plus pos control and NTC)
# from the snounou plates to confirm Pv 18s with sanger sequencing
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



adult_snounou <- readxl::read_excel(path = "/Volumes/share/1. Data/1. Raw Data/Adult_Pv18s/Snounou_Verification_v4.xlsx", sheet = 2)
adult_snounou <- adult_snounou %>% 
  dplyr::mutate(original_plate = stringr::str_split_fixed(SampleName, "_", n=3)[,1],
                original_plate_loc = stringr::str_split_fixed(SampleName, "_", n=3)[,2],
                barcode = stringr::str_split_fixed(SampleName, "_", n=3)[,3]) %>% 
  dplyr::mutate(barcode = tolower(barcode)) # all lower cases, to remove Case Sensitivity


#---------------
# 96 well plate specs
#---------------
plt_rowlet <- sort(rep(LETTERS[1:8], 12))
plt_colnum <- rep(seq(1:12), 8)

#---------------
# Pulli in recode so this is from the actual 16363 smpls 
#---------------
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
pvpos <- dt[dt$pv18s_fctb == "vivpos",]


#---------------------------------------------------------------------------------------------------------
# Randomly select 94 samples for DUFFY actually being used in this study from the 16,363
#---------------------------------------------------------------------------------------------------------
pvsnounoubarcodes <- sample(pvpos$hivrecode_barcode, size = 94, replace = F) %>% 
  tibble(barcode = .)


#---------------
# Find matches in SN
#---------------
sangersnounou <- inner_join(x=sn, y=pvsnounoubarcodes, by = "barcode")
if(nrow(sangersnounou) != 94){
  stop("There was a barcode matching error. See which one it has a one-off base with as you did in the original data generation")
}

#---------------
# write this out for Nick in a way that is easier to work with 
#---------------


sangersnounou <- sangersnounou %>% 
  dplyr::arrange(plate, rowlett, columnnum) %>% 
  dplyr::mutate(sngplt_rowlet = plt_rowlet[1:nrow(.)],
                sngplt_colnum = plt_colnum[1:nrow(.)],
                snounouname = paste0(plate, "_", rowlett, columnnum),
                snounouname = stringr::str_replace(snounouname, "VivID_", "")) 

sangersnounou %>% 
  write.csv(., file = "WetLabWork/platemapsout/confirmatory_Pv18s_Snounou_sanger_smpls_LONG.csv",
            quote = F, row.names = F)

sangersnounou %>% 
  dplyr::select(c("sngplt_rowlet", "sngplt_colnum", "snounouname")) %>% 
  rbind.data.frame(. ,
                   data.frame(
                     sngplt_rowlet = c("H", "H"),
                     sngplt_colnum =  c(11,12),
                     snounouname = c("PC", "NTC"))
  ) %>% 
  tidyr::spread(., key = "sngplt_colnum", value = "snounouname") %>% 
  write.csv(., file = "WetLabWork/platemapsout/confirmatory_Pv18s_Snounou_sanger_smpls_WIDE.csv",
            quote = F, row.names = F)






