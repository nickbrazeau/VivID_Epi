## .................................................................................
## Purpose: Perm backend for icer mode
##
## Notes: 
## .................................................................................

set.seed(48)
remotes::install_github("nickbrazeau/icer")
library(icer)
library(tidyverse)
library(srvyr)
source("R/00-functions_basic.R")

#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)

# need weighed counts
data <- dtsrvy %>% 
  dplyr::mutate(pfmono =  ifelse(pfldh == 1 & pv18s == 0 & po18s == 0, 1, 0),
                pvmono =  ifelse(pfldh == 0 & pv18s == 1 & po18s == 0, 1, 0),
                pfpv =  ifelse(pfldh == 1 & pv18s == 1 & po18s == 0, 1, 0)
  ) %>% 
  dplyr::summarise(
    "pf" = srvyr::survey_total(x=pfmono),
    "pv" = srvyr::survey_total(x=pvmono),
    "pf/pv" = srvyr::survey_total(x=pfpv)
  ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))



data.vec <- unlist(data)
# have to round 
data.vec <- round(data.vec, 0)

#......................
# make icer model
#......................
surv <- new("surveillance")
surv@denominator <- nrow(dt)
surv@cases <- unlist(data.vec)
surv@casenames <- names(data.vec)

#......................
# run icer model 
#......................
# assuming independence
surv.mle <- icer::cooccurence_mle(obj = surv, 
                                  density_func = icer:::independent,
                                  max_moi = 25, poisson = T, plot = F, boot_iter = 25e3) 

#......................
# save out
#......................
dir.create("analyses/08-nonparam_riskfacts/perm_rets/", recursive = T)
saveRDS(surv.mle, "analyses/08-nonparam_riskfacts/perm_rets/icer_model.rds")