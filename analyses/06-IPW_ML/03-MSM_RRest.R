#----------------------------------------------------------------------------------------------------
# Purpose of this script is to calculate IPTW-Weighted Risk Ratios
# for our marginal structural models
#----------------------------------------------------------------------------------------------------


#.....................................
# Covariates that are unconfounded in expectation
#.....................................
# urbanicity and cluster altitude (historical reasons where a cluster is/where people live. Not a causal factor)
# sex (biological chance)
# age (biological process)



#----------------------------------------------------------------------------------------------------
# TABLE 2
# Generalized Propensity Scores
#----------------------------------------------------------------------------------------------------
# note overwriting risk factors to have scaled models now
#.......................
# Pvivax 
#.......................
pvivrskfctr <- dcdr$column_name[dcdr$risk_factor_model == "y"]
pvivrskfctr_models <- data.frame(outcome = rep("pv18s", length(pvivrskfctr)), 
                                 covar = pvivrskfctr, stringsAsFactors=FALSE)

pvivrskfctr_models$glmlogit <- purrr::pmap(pvivrskfctr_models, .f=fitsvyglmlogit)
pvivrskfctr_models$glmlogit_tidy <- purrr::map(pvivrskfctr_models$glmlog,
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
