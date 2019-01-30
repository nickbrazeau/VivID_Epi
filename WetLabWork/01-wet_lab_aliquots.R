#---------------------------------------------------------------------------
# Purpose of this script is to take Snounou plates to long format
# will also randomly select 10% for HRM confirmation by Sanger Sequencing
# will also randomly select ~200 samples for Pv Amplicon Sequencing
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
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
pvpos <- dt[dt$pv18s_fctb == "viv+",]



#---------------------------------------------------------------------------------------------------------
# Randomly select 10% of samples for DUFFY actually being used in this study from the 16,363
#---------------------------------------------------------------------------------------------------------
pvsangbarcodes <- sample(pvpos$hivrecode_barcode, size = 0.1*nrow(pvpos), replace = F) %>% 
  tibble(barcode = .)


#---------------
# Find matches in SN
#---------------
sangerduffy <- inner_join(x=sn, y=pvsangbarcodes, by = "barcode")
if(nrow(sangerduffy) != floor(0.1*nrow(pvpos))){
  stop("There was a barcode matching error. See which one it has a one-off base with as you did in the original data generation")
}

#---------------
# write this out fro Jolly in a way that is easier to work with 
#---------------


sangerduffy <- sangerduffy %>% 
  dplyr::arrange(plate, rowlett, columnnum) %>% 
  dplyr::mutate(sngplt_rowlet = plt_rowlet[1:nrow(.)],
                sngplt_colnum = plt_colnum[1:nrow(.)],
                snounouname = paste0(plate, "_", rowlett, columnnum),
                snounouname = stringr::str_replace(snounouname, "VivID_", "")) 

sangerduffy %>% 
  write.csv(., file = "WetLabWork/platemapsout/confirmatory_duffy_sanger_smpls_LONG.csv",
            quote = F, row.names = F)

sangerduffy %>% 
  dplyr::select(c("sngplt_rowlet", "sngplt_colnum", "snounouname")) %>% 
  rbind.data.frame(. ,
                   data.frame(
                     sngplt_rowlet = c("E", "E", "E"),
                     sngplt_colnum = c(2,3,4),
                     snounouname = c("PC1", "PC2", "NTC"))
                   ) %>% 
  tidyr::spread(., key = "sngplt_colnum", value = "snounouname") %>% 
  write.csv(., file = "WetLabWork/platemapsout/confirmatory_duffy_sanger_smpls_wide.csv",
            quote = F, row.names = F)







#---------------------------------------------------------------------------------------------
# Randomly select 188 samples (-2 for each plate PC and NTC) of samples actually being used in this study from the 16,363
# for the Pv Amplicon
#---------------------------------------------------------------------------------------------
pvampbarcodes <- sample(pvpos$hivrecode_barcode, size = 188, replace = F) %>% 
  tibble(barcode = .)


#---------------
# Find matches in SN
#---------------
pvamp <- inner_join(x=sn, y=pvampbarcodes, by = "barcode")
if(nrow(pvamp) != 188){
  stop("There was a barcode matching error. See which one it has a one-off base with as you did in the original data generation")
}


#---------------
# Wite out for Easier Aliquots for Jolly
#---------------
plt <- c(rep("Plate 1", 94), rep("Plate 2", 94)) 
pvamp <- pvamp %>% 
  dplyr::arrange(plate, rowlett, columnnum) %>% 
  dplyr::mutate(pvampplt_rowlet = c(plt_rowlet[1:94], plt_rowlet[1:94]),
                pvampplt_colnum = c(plt_colnum[1:94], plt_colnum[1:94]),
                snounouname = paste0(plate, "_", rowlett, columnnum),
                snounouname = stringr::str_replace(snounouname, "VivID_", ""),
                newplate = plt) 

pvamp %>% 
  write.csv(., file = "WetLabWork/platemapsout/Pvamplicon_smpls_LONG.csv",
            quote = F, row.names = F)



pvamp %>% 
  dplyr::filter(newplate == "Plate 1") %>% 
  dplyr::select(c("pvampplt_rowlet", "pvampplt_colnum", "snounouname")) %>% 
  rbind.data.frame(. ,
                   data.frame(
                     pvampplt_rowlet = c("H", "H"),
                     pvampplt_colnum =  c(11,12),
                     snounouname = c("PC", "NTC"))
                   ) %>% 
  tidyr::spread(., key = "pvampplt_colnum", value = "snounouname") %>% 
  write.csv(., file = "WetLabWork/platemapsout/Pvamplicon_smpls_wide_PLATE1.csv",
            quote = F, row.names = F)


pvamp %>% 
  dplyr::filter(newplate == "Plate 2") %>% 
  dplyr::select(c("pvampplt_rowlet", "pvampplt_colnum", "snounouname")) %>% 
  rbind.data.frame(. ,
                   data.frame(
                     pvampplt_rowlet = c("H", "H"),
                     pvampplt_colnum =  c(11,12),
                     snounouname = c("PC", "NTC"))
  ) %>% 
  tidyr::spread(., key = "pvampplt_colnum", value = "snounouname") %>% 
  write.csv(., file = "WetLabWork/platemapsout/Pvamplicon_smpls_wide_PLATE2.csv",
            quote = F, row.names = F)





