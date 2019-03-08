#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")
dt <- dt %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum))

#......................
# dependencies
#......................
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
library(tidyverse)
library(broom)
library(survey)

#----------------------------------------------------------------------------------------------------
# Sanity Check 1 -- Weights from AR and PR are not the same
#----------------------------------------------------------------------------------------------------

# per the DHS website
# https://dhsprogram.com/data/Using-DataSets-for-Analysis.cfm#CP_JUMP_14042
# like other variables in DHS datasets, decimal points are not included in the weight variable. Analysts need to divide the sampling weight they are using by 1,000,000
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6,
                hiv05_wi = hiv05/1e6) 
weightsdf <- dt %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n=n(),
                   prweightmean = mean(hv005_wi),
                   arweightmean = mean(hiv05_wi))
head(weightsdf)


#----------------------------------------------------------------------------------------------------
# Sanity Check 2 -- DHS Sampling design and strata
#----------------------------------------------------------------------------------------------------
# https://userforum.dhsprogram.com/index.php?t=msg&goto=179&S=Google

sum(dt$hv021 == dt$hv001) == nrow(dt)

sum(!duplicated(dt$hv023))
sum(!duplicated(dt$adm1name))

#----------------------------------------------------------------------------------------------------
# Sanity Check 3 -- Weights from Survey Package versus glm versus GEE
#----------------------------------------------------------------------------------------------------
# get recoded data now
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")

options(survey.lonely.psu="certainty")
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, strata = hv023, weights = hv005_wi)
dtsrvy_nostrat <- dt %>% srvyr::as_survey_design(ids = hv001, weights = hv005_wi)


s <- survey::svyglm(pv18s ~ 1, 
                    design = dtsrvy,
                    family = binomial(link="logit"))

t <- survey::svyglm(pv18s ~ 1, 
                    design = dtsrvy_nostrat,
                    family = binomial(link="logit"))

u <-    glm(pv18s ~ 1, 
            data = dt,
            weights = hv005_wi,
            family = binomial(link="logit"))

v <-    geepack::geeglm(pv18s ~ 1, 
                        data = dt,
                        weights = hv005_wi,
                        id = hv023,
                        family = binomial(link="logit"))


broom::tidy(s)
broom::tidy(t)
broom::tidy(u)
broom::tidy(v)









