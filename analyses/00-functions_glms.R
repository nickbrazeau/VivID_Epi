
#----------------------------------------------------------------------------------------------------
# glms
#----------------------------------------------------------------------------------------------------

fitglm <- fit_model <- function(outcome, covar){
  
  eq <- as.formula(paste0(outcome, "~", covar))
  ret <- glm(eq,
             data = dt,
             family=binomial(link="logit"))
  
  return(ret)
  
}

fitglm_prov <- fit_model <- function(outcome, covar, data){
  
  eq <- as.formula(paste0(outcome, "~", covar, "+ (1|adm1name)"))
  ret <- glm(eq,
             data = dt,
             family=binomial(link="logit"))
  
  return(ret)
  
}