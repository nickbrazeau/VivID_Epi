#----------------------------------------------------------------------------------------------------
# Purpose of this script is to do basic bivariate analyses
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/R/00-functions_glms.R") 
library(tidyverse)
library(survey)
library(srvyr) #wrap the survey package in dplyr syntax
devtools::install_github("kaz-yos/tableone")
library(tableone)
library(stargazer)
library(nlme)

options("survey.lonely.psu"="certainty")


#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
options(survey.lonely.psu="certainty")
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, strata = hv023, weights = hv005_wi)

#......................
# identify covariates
#......................
vars <- colnames(dt)[grepl("_fctm|_fctb|_scaled", colnames(dt))]
vars <- vars[vars != "pv18s_fctb_ind"]

#......................
# bivariate analyses, covar vs. case
#......................
pvtbl1 <- tableone::svyCreateTableOne(data=dtsrvy)






#----------------------------------------------------------------------------------------------------
# Parametric, Bivariate Analysis
# Odds Ratios with Pv as the outcome
# note, that svyglm is really doing GEE
#----------------------------------------------------------------------------------------------------
model_parameters <- data.frame(outcome = rep("pv18s", length(vars)), 
                               covar = vars, stringsAsFactors=FALSE)

model_parameters$glmlogit <- purrr::pmap(model_parameters, .f=fitsvyglm)
model_parameters$glmlogit_tidy <- purrr::map(model_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))
options(scipen=999)
est <- model_parameters$glmlogit_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2)


m1 <- survey::svyglm(pv18s ~ hv025_fctb_clst, 
                           design = dtsrvy,
                           family = quasibinomial(link = "log"))
broom::tidy(m1, exponentiate = T, conf.int=TRUE)

contrast::contrast(m1,
         a=list(hv025_fctb_clst = "rural"),
         b=list(hv025_fctb_clst = "urban"))









pfmodel_parameters <- data.frame(outcome = rep("pfldh", length(vars)), 
                               covar = vars, stringsAsFactors=FALSE)

pfmodel_parameters$glmlogit <- purrr::pmap(pfmodel_parameters, .f=fitsvyglm)
pfmodel_parameters$glmlogit_tidy <- purrr::map(pfmodel_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))


pomodel_parameters <- data.frame(outcome = rep("po18s", length(vars)), 
                                 covar = vars, stringsAsFactors=FALSE)

pomodel_parameters$glmlogit <- purrr::pmap(pomodel_parameters, .f=
                                             )
pomodel_parameters$glmlogit_tidy <- purrr::map(pomodel_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))
poest <- pomodel_parameters$glmlogit_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2)


m1 <- survey::svyglm(pv18s ~ hml20_fctb_ind + pfldh_fctb_ind, 
                     design = dtsrvy, 
                     family = quasibinomial("log"))
broom::tidy(m1, exponentiate = T, conf.int = T)

m2 <- survey::svyglm(pv18s ~ hml20_fctb_ind + pfldh_fctb_ind + hml20_fctb_ind*pfldh_fctb_ind, 
                     design = dtsrvy, 
                     family = quasibinomial("log"))
broom::tidy(m2, exponentiate = T, conf.int = T)

#----------------------------------------------------------------------------------------------------
# Playgroun
#----------------------------------------------------------------------------------------------------
#...............
# how strong of an effect is that pv outlier
#..............
# not bad bc so few points in cluster 81 even though high prev of viv

m1 <- glm(pv18s ~ hml20_fctb_ind + pfldh + hml20_fctb_ind*pfldh,
          data = dt,
          family = binomial(link = "logit"))
broom::tidy(m1, exponentiate=TRUE, conf.int=TRUE)
contrast::contrast(m1,
         a=list(hml20_fctb_ind = "yes", pfldh = 1),
         b=list(hml20_fctb_ind = "no", pfldh = 1))


#----------------------------------------------------------------------------------------------------
# Save Objects for Reports
#----------------------------------------------------------------------------------------------------
out <- "~/Documents/GitHub/VivID_Epi/reports/report_obj"
if(!dir.exists(out)){dir.create(out)}

save(pftbl1, pvtbl1, model_parameters, file = paste0(out, "/", "03-Uni_Bivar_analyses.rda"))

