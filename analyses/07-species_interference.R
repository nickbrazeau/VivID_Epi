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
    pf = srvyr::survey_total(x=pfmono),
    pv = srvyr::survey_total(x=pvmono),
    po = srvyr::survey_total(x=pomono),
    pfpv = srvyr::survey_total(x=pfpv),
    pfpo = srvyr::survey_total(x=pfpo),
    pvpo = srvyr::survey_total(x=pvpo),
    pfpvpo = srvyr::survey_total(x=pfpvpo)
  ) %>% 
  dplyr::select(-c(dplyr::ends_with("_se")))


 
data.vec <- unlist(data)
# have to round 
real <- data.frame("variable"=c("pf/po","pf/pv","pf/po/pv","pf","po","po/pv","pv"),
                  "value"=c(74,150,1,4667,47,1,326))

# assuming independence
ret.ind <- icer::cooccurence_test(real,
                              boot_iter = 1e4)
 

# assuming interference
ret.interference <- icer::cooccurence_test(real,
                                           density_func = icer:::interference, 
                                           k_12 = 0.5, k_13 = 2, k_23 = 1,
                                           boot_iter = 1e4)


svglite::svglite(file = "~/Desktop/interference.svg", width = 11, height = 16)
cowplot::plot_grid(ret.ind$plot$plot, 
                   ret.interference$plot$plot, 
                   ncol = 1)

graphics.off()










