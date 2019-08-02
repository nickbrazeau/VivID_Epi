#----------------------------------------------------------------------------------------------------
# Purpose of this script is to generate final tables for Publication
#----------------------------------------------------------------------------------------------------
library(tidyverse)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/R/00-functions_epi.R") 


#----------------------------------------------------------------------------------------------------
# MAIN TABLES
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# SUPPLEMENTAL TABLES
#----------------------------------------------------------------------------------------------------
#----------------------------------------------
# Table 1: Case/NonCase Dist
#----------------------------------------------
load("results/bivariate_model_results.rda")

readr::write_csv(x = pvivtbl1df, 
                 path = "results/tables/Tbl1_col1_pvivtbl1.csv",
                 na = "")

readr::write_csv(x = pfaltbl1df, 
                 path = "results/tables/Tbl1_col2_pfaltbl1.csv",
                 na = "")

readr::write_csv(x = casestbl1df, 
                 path = "results/tables/Tbl1_col3_casesbl1.csv",
                 na = "")



#----------------------------------------------
# Supp Table 2: Bivariate, Confounded OR
#----------------------------------------------

readr::write_csv(x = pvivriskfactortable, 
                 path = "results/tables/Supp1_pvivriskfactortable.csv",
                 na = "")

readr::write_csv(x = pfalriskfactortable, 
                 path = "results/tables/Supp2_pfalriskfactortable.csv",
                na = "")









