#----------------------------------------------------------------------------------------------------
# Purpose of this script is to do basic bivariate analyses
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/R/00-functions_epi.R") 
library(tidyverse)
library(survey)
library(srvyr) #wrap the survey package in dplyr syntax
devtools::install_github("kaz-yos/tableone")
library(tableone)

options(scipen=999)

#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "internal_datamap_files/DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))
dtsrvy <- makecd2013survey(survey = dt)
dtnogeo <- dt
sf::st_geometry(dtnogeo) <- NULL
dtsrvy <- makecd2013survey(survey = dtnogeo)


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
# TABLE ONE -- bivariate analyses, descriptive
#----------------------------------------------------------------------------------------------------
#......................
# identify risk factors
#......................
pvivrskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y"]
pfalrskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y" & dcdr$column_name != "pfldh_fctb" & dcdr$column_name != "pfldh_cont_clst" ]
pfalrskfctr <- c("pv18s_fctb", pfalrskfctr)


#.......................
# Pvivax 
#.......................
pvivtbl1 <- tableone::svyCreateTableOne(
                                      data = dtsrvy,
                                      strata = "pv18s_fctb",
                                      vars = pvivrskfctr,
                                      includeNA = T,
                                      test = F)
#.......................
# Pfalciparum 
#.......................
pfaltbl1 <- tableone::svyCreateTableOne(
  data = dtsrvy,
  strata = "pfldh_fctb",
  vars = pfalrskfctr,
  includeNA = T,
  test = F)


#----------------------------------------------------------------------------------------------------
# TABLE 2
# Parametric, Bivariate Analysis
# note, that svyglm is really doing GEE
#----------------------------------------------------------------------------------------------------
#.......................
# Pvivax 
#.......................
pvivrskfctr <- dcdr$column_name[dcdr$risk_factor_model == "y"]
pvivrskfctr_models <- data.frame(outcome = rep("pv18s", length(pvivrskfctr)), 
                               covar = pvivrskfctr, stringsAsFactors=FALSE)

pvivrskfctr_models$glmlogit <- purrr::pmap(pvivrskfctr_models, .f=fitsvyglmlogit)
pvivrskfctr_models$glmlogit_tidy <- purrr::map(pvivrskfctr_models$glmlogit,
                                             .f=function(x){
                                               broom::tidy(x, exponentiate=TRUE, conf.int=TRUE)}
                                             )
pvivrskfctr_est <- pvivrskfctr_models$glmlogit_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2)


#.......................
# Pfalciparum 
#.......................
pfalrskfctr <- dcdr$column_name[dcdr$risk_factor_model == "y" & dcdr$column_name != "pfldh_fctb"]
pfalrskfctr <- c("pv18s_fctb", pfalrskfctr)
pfalrskfctr_models <- data.frame(outcome = rep("pfldh", length(pfalrskfctr)), 
                               covar = pfalrskfctr, stringsAsFactors=FALSE)

pfalrskfctr_models$glmlogit <- purrr::pmap(pfalrskfctr_models, .f=fitsvyglmlogit)
pfalrskfctr_models$glmlogit_tidy <- purrr::map(pfalrskfctr_models$glmlogit,
                                             .f=function(x){
                                               broom::tidy(x, exponentiate=TRUE, conf.int=TRUE)}
)
pfalrskfctr_est <- pfalrskfctr_models$glmlogit_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2)




#----------------------------------------------------------------------------------------------------
# Combine Table 1 and 2
#----------------------------------------------------------------------------------------------------
#.......................
# Pvivax
#.......................
pvivtbl1df <- tableone2dataframe(pvivtbl1, columnnames = c("Covariates",
                                                           "Pvivax-Negative",
                                                           "Pvivax-Positive",
                                                           "matchcol"))
pvivtbl1df <- dcdr %>% 
  dplyr::select(c("column_name", "var_label")) %>% 
  dplyr::rename(matchcol = column_name) %>% 
  dplyr::left_join(pvivtbl1df, ., by = "matchcol") %>% 
  dplyr::select(c("var_label", dplyr::everything()))

# Remember, all continuous variables are scaled in the models but not in the original distributions
pvivtbl1df <- pvivtbl1df %>% 
  dplyr::mutate(matchcol = gsub("_cont_clst", "_cont_scale_clst", matchcol),
                matchcol = gsub("_cont$", "_cont_scale", matchcol))


pvivriskfactortable <- mergetableone2table(tableonedf = pvivtbl1df,
                                           tabletwoestdf = pvivrskfctr_est)


#.......................
# Pfalciparum 
#.......................
pfaltbl1df <- tableone2dataframe(pfaltbl1, columnnames = c("Covariates",
                                                           "Pfalciparum-Negative",
                                                           "Pfalciparum-Positive",
                                                           "matchcol"))

pfaltbl1df <- dcdr %>% 
  dplyr::select(c("column_name", "var_label")) %>% 
  dplyr::rename(matchcol = column_name) %>% 
  dplyr::left_join(pfaltbl1df, ., by = "matchcol") %>% 
  dplyr::select(c("var_label", dplyr::everything()))

# Remember, all continuous variables are scaled in the models but not in the original distributions
pfaltbl1df <- pfaltbl1df %>% 
  dplyr::mutate(matchcol = gsub("_cont_", "_cont_scale_", matchcol),
                matchcol = gsub("_cont$", "_cont_scale", matchcol))


pfalriskfactortable <- mergetableone2table(tableonedf = pfaltbl1df,
                                           tabletwoestdf = pfalrskfctr_est)



#----------------------------------------------------------------------------------------------------
# Save out
#----------------------------------------------------------------------------------------------------
save(pvivtbl1, pfaltbl1, # table one output 
     pvivrskfctr_models, pfalrskfctr_models, # model datatframes
     pvivriskfactortable, pfalriskfactortable, # final out table for report
     file = "results/bivariate_model_results.rda")


#----------------------------------------------------------------------------------------------------
# Playground
# #----------------------------------------------------------------------------------------------------
# 
# 
# 
# m1 <- survey::svyglm(pv18s ~ hv025_fctb_clst, 
#                            design = dtsrvy,
#                            family = quasibinomial(link = "log"))
# broom::tidy(m1, exponentiate = T, conf.int=TRUE)
# 
# contrast::contrast(m1,
#          a=list(hv025_fctb_clst = "rural"),
#          b=list(hv025_fctb_clst = "urban"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# pfmodel_parameters <- data.frame(outcome = rep("pfldh", length(vars)), 
#                                covar = vars, stringsAsFactors=FALSE)
# 
# pfmodel_parameters$glmlogit <- purrr::pmap(pfmodel_parameters, .f=fitsvyglm)
# pfmodel_parameters$glmlogit_tidy <- purrr::map(pfmodel_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))
# 
# 
# pomodel_parameters <- data.frame(outcome = rep("po18s", length(vars)), 
#                                  covar = vars, stringsAsFactors=FALSE)
# 
# pomodel_parameters$glmlogit <- purrr::pmap(pomodel_parameters, .f=
#                                              )
# pomodel_parameters$glmlogit_tidy <- purrr::map(pomodel_parameters$glmlogit, .f=function(x) broom::tidy(x, exponentiate=TRUE, conf.int=TRUE))
# poest <- pomodel_parameters$glmlogit_tidy %>% 
#   bind_rows() %>% filter(term != "(Intercept)") %>% 
#   mutate_if(is.numeric, round, 2)
# 
# 
# m1 <- survey::svyglm(pv18s ~ hml20_fctb_ind + pfldh_fctb_ind, 
#                      design = dtsrvy, 
#                      family = quasibinomial("log"))
# broom::tidy(m1, exponentiate = T, conf.int = T)
# 
# m2 <- survey::svyglm(pv18s ~ hml20_fctb_ind + pfldh_fctb_ind + hml20_fctb_ind*pfldh_fctb_ind, 
#                      design = dtsrvy, 
#                      family = quasibinomial("log"))
# broom::tidy(m2, exponentiate = T, conf.int = T)
# 
# #...............
# # how strong of an effect is that pv cluster 81
# #..............
# 
# 
# #...............
# # look at household level
# #..............
# houseindex <- dt %>% 
#   group_by(houseid) %>% 
#   summarise(n = n(), cases = sum(pv18s), houseinfctperc = cases/n)
# hist(houseindex$houseinfctperc)
# 
# 
# #...............
# # treating bednets as a confounder (if pf is a mediator)
# #..............
# m1 <- survey::svyglm(pv18s ~ pfldh_fctb + hml20_fctb + pfldh_fctb*hml20_fctb, 
#                      design = dtsrvy,
#                      family = quasibinomial("logit"))
# 
# broom::tidy(m1, exponentiate = T, conf.int = T)
# 
# m2 <- glm(pv18s ~ pfldh_fctb + hml20_fctb + pfldh_fctb*hml20_fctb, 
#                      data = dt,
#                      family = quasibinomial("logit"))
# 
# broom::tidy(m2, exponentiate = T, conf.int = T)
# 
# contrast::contrast(m2,
#          a=list(hml20_fctb = "yes", pfldh_fctb = "falneg"),
#          b=list(hml20_fctb = "no", pfldh_fctb = "falneg"))
# 
# 
# 
# #----------------------------------------------------------------------------------------------------
# # Save Objects for Reports
# #----------------------------------------------------------------------------------------------------
# out <- "~/Documents/GitHub/VivID_Epi/reports/report_obj"
# if(!dir.exists(out)){dir.create(out)}
# 
# save(pftbl1, pvtbl1, model_parameters, file = paste0(out, "/", "03-Uni_Bivar_analyses.rda"))
# 
