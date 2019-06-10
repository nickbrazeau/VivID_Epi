#-------------------
# Import data 
#-------------------
# source("analyses/01-data_import_dhs_qpcr.R")
# source("analyses/01-data_import_ecological.R")
# source("analyses/01-data_import_hotosm.R")
# source("analyses/01-data_import_physical_mapaesthetics.R")

#-------------------
# Data Wrangle
#-------------------
source("analyses/02-data_wrangling_actuse.R")
source("analyses/02-data_wrangling_urbanicity.R")
source("analyses/02-data_wrangling_wealth.R")
source("analyses/02-data_wrangling.R")

#-------------------
# Modeling Due Diligence
#-------------------
source("analyses/03-covar_assoc.R")


#-------------------
# Bivariate Models
#-------------------
source("analyses/04-Uni_Bivar_Analyses.R")



#-------------------
# Bring it home
#-------------------
source("999-Generate_Final_Figures.R")
source("999-Generate_Final_Tables.R")















