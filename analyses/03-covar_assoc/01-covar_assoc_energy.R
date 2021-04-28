#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get the associations
# between covariates using slurm. Using this energy
# function as the differing functional forms of covariates is 
# not an issue... let's us compare binary and continuous 
# covariates on the same scale
# 
# https://www.rdocumentation.org/packages/energy/versions/1.7-6/topics/distance%20correlation
# https://cran.r-project.org/web/packages/energy/energy.pdf
# https://stats.stackexchange.com/questions/183572/understanding-distance-correlation-computations
# https://projecteuclid.org/euclid.aoas/1267453933
# http://yunus.hacettepe.edu.tr/~iozkan/eco742/Brownian.html
# https://projecteuclid.org/download/pdfview_1/euclid.aos/1201012979
# https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
# 
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(energy)
set.seed(48, "L'Ecuyer")

#......................
# Import Data
#......................
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map_v3.xlsx", sheet = 1) %>% 
  dplyr::mutate( risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model) )

# grab risk factors
rskfctr <- dcdr %>% 
  dplyr::filter(risk_factor_raw == "y" ) %>% 
  dplyr::select("column_name") %>% 
  unlist(.) %>% 
  unname(.)

# find covars to compare
paramsdf <- t( combn(rskfctr, 2) ) %>% 
  tibble::as_tibble(.) %>% 
  magrittr::set_colnames(c("covar1", "covar2"))



#......................
# Subset to Final Data and remove missingness 
#......................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
dt.ml.cc <- dt  %>% 
  dplyr::select(rskfctr) %>% 
  dplyr::filter(complete.cases(.)) 

# add dt data in for it to find on slurm
paramsdf$data <- lapply(1:nrow(paramsdf), function(x) return(dt.ml.cc))


#......................
# wrapper function
#......................
energy_calc_corr <- function(covar1, covar2, data){
  x1 <- unlist( data[,covar1] )
  x2 <- unlist( data[,covar2] )
  if(is.factor(x1)){
    x1 <- as.numeric(x1)
  }
  if(is.factor(x2)){
    x2 <- as.numeric(x2)
  }
  
  ret <- energy::dcor(x = x1, y = x2)
  
  return(ret)
}


#................................
# run w/ furrr
#................................
future::plan("multicore")
paramsdf <- paramsdf %>% 
  dplyr::mutate(dcor = furrr::future_pmap_dbl(., energy_calc_corr))

# save out 
saveRDS(paramsdf, file = "results/covars_collinearity.RDS")


