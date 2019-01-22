#----------------------------------------------------------------------------------------------------
# Purpose of this script is to do basic bivariate analyses
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions_glms.R") 
library(tidyverse)
library(srvyr) #wrap the survey package in dplyr syntax
devtools::install_github("kaz-yos/tableone")
library(tableone)
library(stargazer)
library(nlme)

#......................
# Import Data
#......................
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")


#----------------------------------------------------------------------------------------------------
# Table One for Pv
#----------------------------------------------------------------------------------------------------
vars <- colnames(dt)[grepl("_fctm|_fctb|_cont", colnames(dt))]
# put drop weights and household id and continous data
vars <- vars[!vars %in% c("hv005_cont", "hiv05_cont", "hhid_fctm", "hvdate_cont", "pv18s_fctb")]
pvtbl1 <- tableone::CreateTableOne(data=dt, 
                                   strata = "pv18s_fctb", 
                                   vars = vars,
                                   includeNA = T)

#----------------------------------------------------------------------------------------------------
# Table One for Pf
#----------------------------------------------------------------------------------------------------
vars <- colnames(dt)[grepl("_fctm|_fctb|_cont", colnames(dt))]
# put drop weights and household id and continous data
vars <- vars[!vars %in% c("hv005_cont", "hiv05_cont", "hhid_fctm", "hvdate_cont", "pfldh_fctb")]
pftbl1 <- tableone::CreateTableOne(data=dt, 
                                   strata = "pfldh_fctb", 
                                   vars = vars,
                                   includeNA = T)


dt %>% 
  group_by(hvyrmnth_fctm) %>% 
  summarise(n=n())



#----------------------------------------------------------------------------------------------------
# Parametric, Bivariate Analysis
# Odds Ratios with Pv as the outcome
#----------------------------------------------------------------------------------------------------
covars <- colnames(dt)[grepl("_fctm|_fctb|_cont", colnames(dt))]

covars <- covars[!covars %in% c("hv005_cont", "hiv05_cont", "hhid_fctm", "hvdate_cont", "pv18s_fctb")]
model_parameters <- data.frame(outcome = rep("pv18s", length(covars)), 
                               covar = covars, stringsAsFactors=FALSE)

model_parameters$glmlogit <- purrr::pmap(model_parameters, .f=fitglm)
model_parameters$glmlogit_tidy <- purrr::map(model_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))


#----------------------------------------------------------------------------------------------------
# Parametric, Bivariate Analysis, Province Level Random Effect
# RANDOM EFFECT Odds Ratios with Pv as the outcome
#----------------------------------------------------------------------------------------------------
dt$adm1name <- as.numeric(factor(dt$adm1name))


model_parameters$glmlogit_mlm <- purrr::pmap(model_parameters[,1:2], .f=fitglm_prov)
model_parameters$glmlogit_mlm_tidy <- purrr::map(model_parameters$glmlogit_mlm, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))




#----------------------------------------------------------------------------------------------------
# Save Objects for Reports
#----------------------------------------------------------------------------------------------------
out <- "~/Documents/GitHub/VivID_Epi/reports/report_obj"
if(!dir.exists(out)){dir.create(out)}

save(pftbl1, pvtbl1, model_parameters, file = paste0(out, "/", "03-Uni_Bivar_analyses.rda"))

