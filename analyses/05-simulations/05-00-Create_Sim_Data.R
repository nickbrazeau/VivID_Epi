#----------------------------------------------------------------------------------------------------
# Purpose of this script is simulate data from the CD2013 framework 
# which is to say that I want to perserve the spatial relationships between clusters
# and have similar Ns for clusters 
#----------------------------------------------------------------------------------------------------
#......................
# Import Data and dependencies
#......................
source("analyses/00-functions.R") 
library(tidyverse)

load("data/vividepi_raw.rda")
dt <- merge_pr_plsmdm_gemtdt(pr = arpr, plsmdm = panplasmpcrres, ge = ge)

#......................
# drop to important covariates for sim
#......................
dt <- dt %>% 
  select(c("hivrecode_barcode", "hv002", "pv18s", "adm1dhs", "adm1name", "dhsregco", "dhsregna",
           "latnum", "longnum", "geometry"))

#..........................................................
# Create simulation function for independent covariates
#..........................................................
 
# TEMPORARY UNTIL LAST SNOUNOU 
dt$pv18s[dt$hivrecode_barcode == "e3a6n"] <- 1
outcome_prev = mean(dt$pv18s)


simulate_cd2013_indcov <- function(outcome_prev = 0.031,
                            exp_prev = list(),
                            n_covar = 10, 
                            n_covar_pred = 2,
                            n_covar_pred_effect = list()){
  
  
  #  I have assumed that the n_covar_pred_effect is p(x=0|Y)/p(x=0) via bayes
  
  
  #......................
  # Error handling
  #......................
  if(n_covar_pred != length(n_covar_pred_effect)){
    stop("You must provide an effect for each of the covariates that you want to predict the outcome.")
  }
  if(n_covar_pred != length(exp_prev)){
    stop("You must provide a prevalence of exposure for each of the covariates that you want to predict the outcome.")
  }
  
  if(any(unlist(n_covar_pred_effect) < 0 | unlist(n_covar_pred_effect) > 1 )){
    stop("The effect of the covariate must be a prob between 0 and 1.")
  }
  
  #......................
  # setup assoc params
  #......................
  dfout <- tibble::tibble(Y = rbinom(n = nrow(dt), 1, prob = outcome_prev))
  expassoc_store <- list()
  
  for(i in 1:length(n_covar_pred_effect)){
    # P(Y|X=1) = P(Y) - P(Y|X=0)
    Pyx1 <- mean(dfout$Y) - mean(dfout$Y)*n_covar_pred_effect[[i]]
    
    # store 
    expassoc <- Pyx1 - mean(dfout$Y)*n_covar_pred_effect[[i]]
    expassoc_store <- unlist(append(expassoc_store, expassoc))
    
    # P(X|Y=1) = P(Y,X)/P(Y) = P(Y|X=1)*P(X=1) / P(Y)
    Pxy1 <- (Pyx1 * exp_prev[[i]])/outcome_prev
    # P(X|Y=0) = P(X) - P(Y|X=1)
    Pxy0 <- exp_prev[[i]] - Pxy1
    
    temp <- dfout %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(col = 
             ifelse(Y == 1, rbinom(n=1, 1, Pxy1), rbinom(n=1, 1, Pxy0))
           ) 
    dfout <- dplyr::bind_cols(dfout, temp[,2])


    
  }
  
  
  
    
  #......................
  # setup random params
  #......................
  iters <- 1e4 # does this make sense to be a global option?
  # generate some random params
  # thanks hadley, https://r4ds.had.co.nz/iteration.html#invoking-different-functionss
  params <- lapply(1:iters, function(x){
    
    ret <- tribble(
      ~f,      ~params,
      "runif", list(min = runif(1, -100, 0), max =  runif(1, 0, 100)),
      "rnorm", list( runif(1, 0,100)),
      "rpois", list(lambda =  runif(1, 0,100)),
      "rbinom", list(size = 1, prob = runif(1, 0,1))
    )
    return(ret)
  }) %>% dplyr::bind_rows(.)
  
  # extract random params
  params <- params[sample(x = 1:iters, size = n_covar - n_covar_pred, replace = F), ]
  
  
  # invoke params
  params_out <- params %>% 
    dplyr::mutate(sim = invoke_map(f, params, n = nrow(dt))) %>% 
    dplyr::select(-c("f")) %>% 
    tidyr::spread(., key="params", value = "sim") %>% 
    tidyr::unnest()
  
  #......................
  # bring it all home
  #......................
  dfout <- dplyr::bind_cols(dfout, params_out) %>% 
    magrittr::set_colnames(., c("Y", paste0("x", 1:n_covar))) 
  
  ret <- list(
    simdata = dfout, 
    exp_assoc = expassoc_store
  )
    
  return(ret)
  
}


#..........................................................
# Create simulation function for hierarchical model/covariates
#..........................................................

# TEMPORARY UNTIL LAST SNOUNOU 
dt$pv18s[dt$hivrecode_barcode == "e3a6n"] <- 1
outcome_prev = mean(dt$pv18s)


simulate_cd2013_hierarchical <- function(outcome_prev = 0.031,
                                   exp_prev = list(),
                                   n_covar = 10, 
                                   n_covar_pred = 2,
                                   n_covar_pred_effect = list(),
                                   cluster_hier = T,
                                   house_hier = T){
  

# I don't know if it is enough to start with just a group_by
# need to think about making these dependent...potentially 
# regressing onto each other 
# or dropping in a network weighted by road distance matrix


}


