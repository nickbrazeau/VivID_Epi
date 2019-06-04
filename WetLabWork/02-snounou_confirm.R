#---------------------------------------------------------------------------
# Purpose of this script is to randomly select 95 samples
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




#---------------
# 96 well plate specs
#---------------
plt_rowlet <- sort(rep(LETTERS[1:8], 12))
plt_colnum <- rep(seq(1:12), 8)




#---------------
# Subset to plates that we still have
#---------------
sn <- sn %>% 
  dplyr::filter(plate %in% c("VivID_M021-M041", "VivID_M042-M062_plate_1", "VivID_M042-M062_plate_2", "VivID_M063-M084", 
                             "VivID_M124-M143", "VivID_M144-M163", "VivID_M164-M183"))


#---------------
# Pull in recode so this is from the actual smpls being used (and know they were snounou positive already)
#---------------
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
pvpos <- dt[dt$pv18s_fctb == "vivpos",] %>% 
  dplyr::rename(barcode = hivrecode_barcode)



# note, this already isn't a truly random sample since we have plates that are 
# missing. So going to ignore the issue with some samples being lost to barcode
# typos as well and just select 95 samples that I know are being used by my study
# and we have sample left to sanger sequene

snounouret <- dplyr::inner_join(pvpos[,"barcode"], sn)
snounouret <- sample(snounouret$barcode, size = 95, replace = F)
sangersnounou <- sn %>% 
  dplyr::filter(barcode %in% snounouret)

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
                     sngplt_rowlet = c("H"),
                     sngplt_colnum =  c(12),
                     snounouname = c("PC"))
  ) %>% 
  tidyr::spread(., key = "sngplt_colnum", value = "snounouname") %>% 
  write.csv(., file = "WetLabWork/platemapsout/confirmatory_Pv18s_Snounou_sanger_smpls_WIDE.csv",
            quote = F, row.names = F)





# now write out barcode for the Sanger Submission
sangersnounou %>% 
  dplyr::select(c("sngplt_rowlet", "sngplt_colnum", "barcode")) %>% 
  rbind.data.frame(. ,
                   data.frame(
                     sngplt_rowlet = c("H"),
                     sngplt_colnum =  c(12),
                     barcode = c("PC"))
  ) %>% 
  tidyr::spread(., key = "sngplt_colnum", value = "barcode") %>% 
  write.csv(., file = "WetLabWork/platemapsout/Sanger_Submission_Sheet_confirmatory_Pv18s_Snounou_sanger_smpls_WIDE.csv",
            quote = F, row.names = F)




