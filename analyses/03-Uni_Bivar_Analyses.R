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
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))
sf::st_geometry(dt) <- NULL
dtsrvy <- makecd2013survey(survey = dt)

#----------------------------------------------------------------------------------------------------
# Basic Descriptive Statistics
#----------------------------------------------------------------------------------------------------
# national prevalence
sumnums <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::summarise(n = srvyr::survey_total(count, vartype = c("se", "ci")),
                   pvn = srvyr::survey_total(pv18s, vartype = c("se", "ci")),
                   pfn = srvyr::survey_total(pfldh, vartype = c("se", "ci")))

drcpv <- survey::svyglm(pv18s ~ 1, design = dtsrvy)
broom::tidy(drcpv, conf.int = T)

# national prevalence
drcpf <- survey::svyglm(pfldh ~ 1, design = dtsrvy)
broom::tidy(drcpf, conf.int = T)

# # cluster-level prevalence, because everyone is weighted the same in the cluster, don't use weights here
# # note, the numerators will be slightly different (e.g. N) but the denomminators adjust for this
# # Going to report whole numbers/unadjusted for clusters
# clst <- dt %>% 
#   dplyr::mutate(count = 1) %>% 
#   dplyr::group_by(hv001) %>% 
#   dplyr::summarise(n = n(), 
#                    pv18sn = sum(pv18s), 
#                    pv18sprev = mean(pv18sn),
#                    pv18sse = sqrt(pv18sn * pv18sprev * (1 - pv18sprev)) / sqrt(pv18sn),
#                    pv18sprevU95 = pv18sprev + 1.96 * pv18sse,
#                    pv18sprevL95 = pv18sprev - 1.96 * pv18sse,
#                    
#                    pfldhn = sum(pfldh), 
#                    pfldhprev = mean(pfldhn),
#                    pfldhse = sqrt(pfldhn * pfldhprev * (1 - pfldhprev)) / sqrt(pfldhn),
#                    pfldhprevU95 = pfldhprev + 1.96 * pfldhse,
#                    pfldhprevL95 = pfldhprev - 1.96 * pfldhse
#                   
#                    )
#                    
# 
# 
# summary(clst$pv18sn)
# summary(clst$pv18sprev)
# summary(clst$pfldhn)
# summary(clst$pfldhprev)


# ignore 95% CI but want the weighted numerator and denominator for the paper 
clstprev <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = srvyr::survey_total(count, vartype = c("se", "ci")),
                   pvn = srvyr::survey_total(pv18s, vartype = c("se", "ci")),
                   pfn = srvyr::survey_total(pfldh, vartype = c("se", "ci")),
                   
                   pvprev = srvyr::survey_mean(pv18s, vartype = c("se", "ci")),
                   pfprev = srvyr::survey_mean(pfldh, vartype = c("se", "ci"))
                   )



# pfldh coinfection
coinfx <- dtsrvy %>% 
  dplyr::mutate(pvpf = ifelse(pv18s == 1 & pfldh == 1, 1, 0)) %>% 
  srvyr::summarise(
    pfpvcoinfx = srvyr::survey_total(pvpf, vartype = "ci")
  )

xtabs(~pv18s + pfldh, data = dt)




#----------------------------------------------------------------------------------------------------
# TABLE ONE -- bivariate analyses, descriptive
#----------------------------------------------------------------------------------------------------
#......................
# identify risk factors
#......................
pvivrskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y"]
pfalrskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y" & dcdr$column_name != "pfldh_fctb" ]
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
# note overwriting risk factors to have scaled models now
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
                                             broom::tidy(x, exponentiate=TRUE, conf.int=TRUE)})


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


pfalriskfactortable <- mergetableone2table(tableonedf = pfaltbl1df,
                                           tabletwoestdf = pfalrskfctr_est)




#----------------------------------------------------------------------------------------------------
# Make Pv Cases, Pf Cases, and Non-Cases Table
#----------------------------------------------------------------------------------------------------
dt.cases <- dt %>% 
  dplyr::mutate(
    case = ifelse(( pv18s == 1 | pfldh ==1 ), 1, 0),
    case_fctb = factor(case, levels = c(0,1), labels = c("noncase", "case"))
  )

dt.cases.srvy <- makecd2013survey(dt.cases)
casesrskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y"]
casestbl1 <- tableone::svyCreateTableOne(
  data = dt.cases.srvy,
  strata = "case_fctb",
  vars = casesrskfctr,
  includeNA = T,
  test = F)


casestbl1df <- tableone2dataframe(casestbl1, columnnames = c("Covariates",
                                                           "Case-Negative",
                                                           "Case-Positive",
                                                           "matchcol"))


#----------------------------------------------------------------------------------------------------
# Save out
#----------------------------------------------------------------------------------------------------
save(pvivtbl1df, pfaltbl1df, # table one output 
     pvivrskfctr_models, pfalrskfctr_models, # model datatframes
     pvivriskfactortable, pfalriskfactortable, # final out table for report
     casestbl1df, # for prettier table 1
     file = "results/bivariate_model_results.rda")

#----------------------------------------------------------------------------------------------------
# playground
#----------------------------------------------------------------------------------------------------
# m1 <- glm(pv18s ~ precip_ann_cont_scale_clst,
#           family = quasibinomial(link="log"),
#           data = dt)
# 
# m1 <- svyglm(pv18s ~ precip_ann_cont_scale_clst,
#           family = binomial(link="logit"),
#           design = dtsrvy)
# broom::tidy(m1, conf.int = T, exponentiate =T)
# broom::tidy(m1)
# confint_tidy(m1, exponentiate=T)





