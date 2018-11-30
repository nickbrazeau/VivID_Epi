#......................
# Import Data
#......................
source("analyses/00-functions.R") 
ge <- readRDS("datasets/CDGE61FL.rds")
load("data/vividepi_raw.rda")
dt <- merge_pr_plsmdm_gemtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge)
#......................
# dependencies
#......................
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
  dplyr::mutate(hv005 = hv005/1e5,
                hiv05 = hiv05/1e5) 
weightsdf <- dt %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n=n(),
                   prweightmean = mean(hv005),
                   arweightmean = mean(hiv05))
head(weightsdf)

#----------------------------------------------------------------------------------------------------
# Sanity Check 2 -- Weights from Survey Package and glm are the same
#----------------------------------------------------------------------------------------------------

#.......................
# survey package approach
#.......................
options(survey.lonely.psu="adjust")
# clusters are weighted (each individual has same weight in cluster)
design <- survey::svydesign( ids= ~hv001, weights= ~hv005, data=dt)
svyglm(pfldh~1,design)
confint(svyglm(pfldh~1,design))


#.......................
# base glm approach
#.......................
mod1 <- glm(pfldh ~ 1, 
            family = binomial(link = "identity"), 
            data=dt, weights = hv005)

broom::tidy(mod1, conf.int=TRUE)


mod2 <- glm(pfldh ~ factor(adm1name), 
            family = binomial(link = "identity"), 
            data=dt, weights = hv005)

broom::tidy(mod2, conf.int=TRUE)

# should compare this with both nest and with contrast statements







