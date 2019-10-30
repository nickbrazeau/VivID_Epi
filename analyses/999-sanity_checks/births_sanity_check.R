births <- readRDS(file = "~/Documents/GitHub/18Fall_EPID799C_RforEpi/data/births_afterhw2.RDS")

#........................................................
# Births Code from R Class (Epi799C that we taught -- on my github)
#........................................................
df <- births %>% 
  dplyr::select(x, wksgest, preterm, preterm_f, pnc5, pnc5_f, mage, raceeth_f, smoker_f) %>%
  dplyr::filter(!is.na(preterm_f) & !is.na(pnc5_f) & !is.na(mage) & !is.na(raceeth_f) & !is.na(smoker_f))
tableone::CreateTableOne(vars = c("mage","raceeth_f","smoker_f"), 
                         strata = "pnc5_f", test = F,
                         data = df)

df$mage2 <- df$mage^2
covars = c("mage", "mage2", "raceeth_f", "smoker_f")
covars.form = paste(covars, collapse = "+")
summary(ps_model <- glm(data = df, pnc5_f ~ mage+mage2+raceeth_f+smoker_f, family=binomial("logit")))
df$ps <- predict(ps_model, type = "response")

#Distribution of propensity score...
summary(df$ps[df$pnc5==1]) #...among mothers with observed early prenatal care
summary(df$ps[df$pnc5==0]) #...among mothers with no observed early prenatal care

#Use propensity score to calculate unstabilized IPTW
df$iptw_u <- ifelse(df$pnc5==1, 1/df$ps, 1/(1-df$ps)) 

#Distribution of weights in the exposed
summary(df$iptw_u[df$pnc5==1]) 
summary(df$iptw_u[df$pnc5==0])
summary(df$iptw_u)


summary(df$iptw_u)
#Distribution of weights in the exposed
p_exposure <- sum(df$pnc5) / nrow(df)
df$iptw_s <- ifelse(df$pnc5==1, p_exposure/df$ps, (1-p_exposure)/(1-df$ps))
summary(df$iptw_s) 


#........................................................
# mycode
#........................................................
library(mlr)
source("R/00-IPTW_functions.R")

task <- makeClassifTask(data = df[,c(covars, "pnc5_f")], target = "pnc5_f")
lrn <- makeLearner("classif.logreg", predict.type = "prob")
lrn <- list(lrn)

nsets <- 1:nrow(mlr::getTaskData(task))
nsmp <- sample(c(nsets), size = length(nsets)/2)
nsmp <- ifelse(nsets %in% nsmp, 1, 2)
valset.list <- split(nsets, factor(nsmp))


out <- mlrwrapSL::SL_crossval_risk_pred(learnerlib = lrn, task = task, valset.list = valset.list)
summary( get_iptw_prob(task = task, ELpreds = out$EL.predictions) )


