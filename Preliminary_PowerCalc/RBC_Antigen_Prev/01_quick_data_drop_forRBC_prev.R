# tidy
library(tidyverse)
# spatial
library(sf)

##----------------------------------##
##            Read-In Data          ## 
##----------------------------------##
dhs13 <- read_csv(file="/Volumes/share/2. Projects/1. Current/AdultDHS_EpidemiologicalAnalysis/alladults_v3_clean.csv", col_names = T)
colnames(dhs13) <- tolower(colnames(dhs13))

##----------------------------------##
##         Sampling Weights         ## 
##----------------------------------##
# per the DHS website
# https://dhsprogram.com/data/Using-DataSets-for-Analysis.cfm#CP_JUMP_14042
# like other variables in DHS datasets, decimal points are not included in the weight variable. Analysts need to divide the sampling weight they are using by 1,000,000
dhs13$sampleweight <- dhs13$hv005/1000000


##----------------------------------##
##           Shape Objects          ## 
##----------------------------------##
DRCprovfiles <- dir("~/Google Drive/PhD_Work/Maps/DRC_GADM/sf_gadm/", full.names = T)
DRC <- lapply(DRCprovfiles, function(x)readRDS(x))
DRCprov <- DRC[[2]]

provdict <- read_tsv(file = "~/Google Drive/PhD_Work/Maps/DRC_GADM/dhs_shnprovin_mapfile_variabledocumentation.tab.txt")

#---------------------------
######################################
###     Cedars Pvivax Results      ###
######################################
#f <- file.choose()
f <- "~/Desktop/Database_Pv_PoScreen_20180629.xlsx"
pv <- readxl::read_excel(path=f, sheet = 1)
colnames(pv) <- tolower(colnames(pv))
# pv[pv == "Undetermined"] # this doesn't work bc of tibble?
pv$pv_pos <- ifelse(pv$`p. vivax 18s (ct)` != "Undetermined", 1, ifelse(is.na(pv$`p. vivax 18s (ct)`), NA, 0))





######################################
###   Tidy Cedars Pvivax Results   ###
######################################
# drop standards
pv <- pv[pv$task == "UNKNOWN",] # subset to non-strandard/ntc
pv <- pv[pv$id != "0",] # missing barcodes
table(pv$pv_pos, useNA = "always")

# make dhs13 comparable plate barcode
dhs13$dnaplate_short <- substr(x=as.character(dhs13$dnaplate), start=nchar(as.character(dhs13$dnaplate))-3, stop=nchar(as.character(dhs13$dnaplate)))

pv$`original plate m#`[pv$`original plate m#` == "M078"] <- "C078"
dhs13$dnaplate_short <- gsub("-", "UNC", dhs13$dnaplate_short)

# sanity check -- any of Cedar's dna plate names missing from dhs
if(FALSE %in% (pv$`original plate m#` %in% dhs13$dnaplate_short)){
  stop("Plate names from Cedar's file don't match DHS")
}


## rename cedar's file to match
pv <- pv[, colnames(pv) %in% c("original plate m#", "id", "pv_pos", "p. vivax 18s (ct)")]
colnames(pv)[1:2] <- c("dnaplate_short", "barcode")
pv$barcode <- toupper(pv$barcode)


m <- inner_join(dhs13, pv, by=c("dnaplate_short", "barcode"))
sum(m$pv_pos, na.rm = T)
sum(pv$pv_pos, na.rm = T)

### this is an issue with the barcode for the plate maps versus the barcodes from the DHS -- because 50 missing

prelim <- left_join(dhs13, pv, by=c("dnaplate_short", "barcode"))
#-------------------------------------------------------------------------







##----------------------------------##
##           Final Objects          ## 
##----------------------------------##
dhs13_short <- prelim %>% 
  dplyr::select(barcode, cluster, shnprovin, pv_pos, sampleweight, drcshp13.latnum, drcshp13.longnum)


save(dhs13_short, DRCprov, provdict, prelim, file = "~/Documents/GitHub/VivID_PhD_Thesis/Preliminary_PowerCalc/RBC_Antigen_Prev/data.RData")


  