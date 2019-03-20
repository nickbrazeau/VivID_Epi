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
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, 
                                         strata = hv023, 
                                         weights = hv005_wi)



#----------------------------------------------------------------------------------------------------
# Basic Descriptive Statistics
#----------------------------------------------------------------------------------------------------
# national prevalence
m1 <- survey::svyglm(pv18s ~ 1, design = dtsrvy)
broom::tidy(m1, conf.int = T)

# cluster-level prevalence
clst <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = srvyr::survey_total(count), 
                 pv18sn = srvyr::survey_total(pv18s), 
                 pv18sprev = srvyr::survey_mean(pv18s, vartype = c("ci"), level = 0.95))

summary(clst$pv18sn)
summary(clst$pv18sprev)


# pfldh coinfection
xtabs(~pv18s + pfldh, data = dt)




#----------------------------------------------------------------------------------------------------
# TABLE ONE
#----------------------------------------------------------------------------------------------------
#......................
# identify covariates
#......................
vars <- colnames(dt)[grepl("_fctm|_fctb|_scaled", colnames(dt))]
vars <- vars[vars %in% c("pv18s_fctb", "pv18s_fctb_sens")]

#......................
# bivariate analyses, covar vs. case
#......................
pvtbl1 <- tableone::svyCreateTableOne(data=dtsrvy)






#----------------------------------------------------------------------------------------------------
# TABLE 2
# Parametric, Bivariate Analysis
# Odds Ratios with Pv as the outcome
# note, that svyglm is really doing GEE
#----------------------------------------------------------------------------------------------------
model_parameters <- data.frame(outcome = rep("pv18s_sens", length(vars)), 
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
# Playground
#----------------------------------------------------------------------------------------------------
#...............
# how strong of an effect is that pv cluster 81
#..............


#...............
# look at household level
#..............
houseindex <- dt %>% 
  group_by(houseid) %>% 
  summarise(n = n(), cases = sum(pv18s), houseinfctperc = cases/n)
hist(houseindex$houseinfctperc)


#...............
# treating bednets as a confounder (if pf is a mediator)
#..............
m1 <- survey::svyglm(pv18s ~ pfldh_fctb + hml20_fctb + pfldh_fctb*hml20_fctb, 
                     design = dtsrvy,
                     family = quasibinomial("logit"))

broom::tidy(m1, exponentiate = T, conf.int = T)

m2 <- glm(pv18s ~ pfldh_fctb + hml20_fctb + pfldh_fctb*hml20_fctb, 
                     data = dt,
                     family = quasibinomial("logit"))

broom::tidy(m2, exponentiate = T, conf.int = T)

contrast::contrast(m2,
         a=list(hml20_fctb = "yes", pfldh_fctb = "falneg"),
         b=list(hml20_fctb = "no", pfldh_fctb = "falneg"))



#----------------------------------------------------------------------------------------------------
# Save Objects for Reports
#----------------------------------------------------------------------------------------------------
out <- "~/Documents/GitHub/VivID_Epi/reports/report_obj"
if(!dir.exists(out)){dir.create(out)}

save(pftbl1, pvtbl1, model_parameters, file = paste0(out, "/", "03-Uni_Bivar_analyses.rda"))

