#----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate if there is interference between
# Pfal and Pviv infections
#----------------------------------------------------------------------------------------------------
remotes::install_github("OJWatson/icer")
library(icer)
source("R/00-functions_basic.R")

#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
dtsrvy <- makecd2013survey(survey = dt)

# need weighed counts
data <- dtsrvy %>% 
  dplyr::mutate(pfmono =  ifelse(pfldh == 1 & pv18s == 0 & po18s == 0, 1, 0),
                pvmono =  ifelse(pfldh == 0 & pv18s == 1 & po18s == 0, 1, 0),
                pomono =  ifelse(pfldh == 0 & pv18s == 0 & po18s == 1, 1, 0),
                pfpv =  ifelse(pfldh == 1 & pv18s == 1 & po18s == 0, 1, 0),
                pfpo =  ifelse(pfldh == 1 & pv18s == 0 & po18s == 1, 1, 0),
                pvpo =  ifelse(pfldh == 0 & pv18s == 1 & po18s == 1, 1, 0),
                pfpvpo =  ifelse(pfldh == 1 & pv18s == 1 & po18s == 1, 1, 0)
  ) %>% 
  dplyr::summarise(
    "pf" = srvyr::survey_total(x=pfmono),
    "pv" = srvyr::survey_total(x=pvmono),
    "po" = srvyr::survey_total(x=pomono),
    "pf/pv" = srvyr::survey_total(x=pfpv),
    "pf/po" = srvyr::survey_total(x=pfpo),
    "pv/po" = srvyr::survey_total(x=pvpo),
    "pf/pv/po" = srvyr::survey_total(x=pfpvpo)
  ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))



data.vec <- unlist(data)
# have to round 
data.vec <- round(data.vec, 0)

# assuming independence
ret.ind <- icer::cooccurence_test(data.vec,
                                  boot_iter = 5e4)


# assuming interference
ret.interference <- icer::cooccurence_test(data.vec,
                                           density_func = icer:::interference, 
                                           k_12 = 0.5, k_13 = 2, k_23 = 1, # arb starting params
                                           boot_iter = 5e4)


#----------------------------------------------------------------------------------------------------
# Save out
#----------------------------------------------------------------------------------------------------
save(ret.ind, ret.interference,
     file = "results/icer_interference_models.rda")





