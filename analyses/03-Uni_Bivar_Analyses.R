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






dtsrvy <- survey::svydesign(ids = ~hv021 + hv002, weights = ~hiv05_wi, data = dt)

s <- survey::svyglm(pv18s ~ hv025_fctb_clst, 
                    design = dtsrvy,
                    family = quasibinomial(link="logit"))

broom::tidy(s, exponentiate = T, conf.int = T)




dt$strata <- paste0(dt$adm1name, dt$hv001)
svy <- dt %>%  srvyr::as_survey_design(strata=strata, weights=hiv05_wi)

u <- survey::svyglm(pv18s ~ hv025_fctb_clst, 
                    design = svy,
                    family = quasibinomial(link="logit"))

broom::tidy(u, exponentiate = T, conf.int = T)




#----------------------------------------------------------------------------------------------------
# Table One for Pv
#----------------------------------------------------------------------------------------------------
vars <- colnames(dt)[grepl("_fctm|_fctb|_scaled", colnames(dt))]
vars <- vars[vars != "pv18s_fctb_ind"]
  
pvtbl1 <- tableone::svyCreateTableOne(data=dtsrvy)

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

covars <- covars[!covars %in% c("hv005_cont", "hiv05_cont", "hhid_fctm", "hvdate_cont", "pv18s_fctb", "pv18sct_cont"
                                )]



model_parameters <- data.frame(outcome = rep("pv18s", length(vars)), 
                               covar = vars, stringsAsFactors=FALSE)

model_parameters$glmlogit <- purrr::pmap(model_parameters, .f=fitsvyglm)
model_parameters$glmlogit_tidy <- purrr::map(model_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))

dt$id <- paste0(dt$hv021, dt$hv002)
m1 <- geepack::geeglm(pfldh ~ hv025_fctb_clst, 
          data = dt, 
          family=binomial(link="logit"),
          id = id
          )

broom::tidy(m1, exponentiate=TRUE, conf.int=TRUE)

#----------------------------------------------------------------------------------------------------
# Parametric, Bivariate Analysis, Province Level Random Effect
# RANDOM EFFECT Odds Ratios with Pv as the outcome
#----------------------------------------------------------------------------------------------------
dt$adm1name <- as.numeric(factor(dt$adm1name))


model_parameters$glmlogit_mlm <- purrr::pmap(model_parameters[,1:2], .f=fitglm_prov)
model_parameters$glmlogit_mlm_tidy <- purrr::map(model_parameters$glmlogit_mlm, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))





#----------------------------------------------------------------------------------------------------
# Playgroun
#----------------------------------------------------------------------------------------------------
#...............
# how strong of an effect is that pv outlier
#..............
# not bad bc so few points in cluster 81 even though high prev of viv

m1 <- glm(pv18s ~ hv025_fctb,
         data = dt,
         family = binomial(link = "logit"))
broom::tidy(t, exponentiate=TRUE, conf.int=TRUE)

m2 <- glm(pv18s ~ hv025_fctb,
         data = dt[dt$hv001 != 81, ],
         family = binomial(link = "logit"))
broom::tidy(u, exponentiate=TRUE, conf.int=TRUE)




#----------------------------------------------------------------------------------------------------
# Save Objects for Reports
#----------------------------------------------------------------------------------------------------
out <- "~/Documents/GitHub/VivID_Epi/reports/report_obj"
if(!dir.exists(out)){dir.create(out)}

save(pftbl1, pvtbl1, model_parameters, file = paste0(out, "/", "03-Uni_Bivar_analyses.rda"))

