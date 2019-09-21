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


#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode.rds")
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

# add dt data in for it to find on slurm
paramsdf$data <- lapply(1:nrow(paramsdf), function(x) return(dt))

# slurm function
slurm_calc_corr <- function(covar1, covar2, data){
  x1 <- data[,covar1]
  x2 <- data[,covar2]
  if(is.factor(x1)){
    x1 <- as.numeric(x1)
  }
  if(is.factor(x2)){
    x2 <- as.numeric(x2)
  }
  
  if(is.factor(x[,2])){
    x[,2] <- as.numeric(x[,2])
  }
  
  ret <- energy::dcor(x = x1, y = x2)
  
  return(ret)
}


#................................
# send it out with rSLURM
#................................


# for slurm on LL
setwd("analyses/05-covar_assoc/")
ntry <- nrow(paramsdf)
sjob <- rslurm::slurm_apply(f = slurm_calc_corr, 
                            params = paramsdf, 
                            jobname = 'covar_corr',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 128),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5-00:00:00"))




