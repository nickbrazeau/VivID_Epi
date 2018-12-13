#----------------------------------------------------------------------------------------------------
# Purpose of this script is to clean and recode the CD2013 PR recode data and potentially the GC data
#----------------------------------------------------------------------------------------------------
source("analyses/00-functions.R") 
library(tidyverse)
devtools::install_github("OJWatson/rdhs", ref="master")
library(rdhs)
#......................
# Import Data
#......................
load("data/vividepi_raw.rda")
dt <- merge_pr_plsmdm_gemtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge)



#......................
# Pulling map file
#......................
# cd2013 was under phase 6
# https://dhsprogram.com/publications/publication-DHSG4-DHS-Questionnaires-and-Manuals.cfm
# recode map https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
mtdt <- as.tibble(rdhs::get_variable_labels(dt))
readr::write_csv(mtdt, path = "internal_datamap_files/covar_names_labels.csv")



# dhs strings that need to be NA
toMatch <- c("missing", "don't know", "don't know brand", "not applicable",
             "inconsistent", "flagged cases")
# HB61, HB62, HB63 have several responses for test not performed -- var won't be included in final analysis

# start
tibble(label = attr(dt$hml1, "label", exact = TRUE),
       response = names(attr(dt$hml1, "labels", exact = TRUE)),
       value = unname(attr(dt$hml1, "labels", exact = TRUE))
      ) %>% 
  dplyr::mutate(recode = dplyr::if_else(grepl(pattern= paste(toMatch,collapse="|"), tolower(value)),
                                              NA, value, "error")
                )

test <- c("missing", "don't know", "don't know brand", "good", "not applicable")
grepl(pattern = paste(toMatch,collapse="|"), test)

apply(dt[,grepl("sh2", colnames(dt))], 2, function(x) sum(!is.na(x))
) 
apply(dt[,grepl("sh23", colnames(dt))], 2, function(x) sum(!is.na(x))
) 
apply(dt[,grepl("hc", colnames(dt))], 2, function(x) sum(!is.na(x))
) 

sum(!is.na(dt$sh240))

t <- if_else(!is.na(dt$ha40), dt$ha40, dt$hb40)


