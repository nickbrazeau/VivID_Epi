#----------------------------------------------------------------------------------------------------
# Purpose of this script is to examine the associations
# between covariates and get a baseline understanding of their joint distributions
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(energy)


#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate( risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model) )

# grab risk factors
rskfctr <- dcdr %>% 
  dplyr::filter(risk_factor_model == "y" ) %>% 
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





