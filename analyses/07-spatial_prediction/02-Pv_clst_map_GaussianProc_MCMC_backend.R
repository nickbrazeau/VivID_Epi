#----------------------------------------------------------------------------------------------------
# Purpose of this script is to produce spatial models
# from prevmap for Pv 
#----------------------------------------------------------------------------------------------------
source("R/00-functions_basic.R")
source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) 
library(raster)
library(PrevMap)
library(drake)
set.seed(48, "L'Ecuyer")
tol <- 1e-3
#......................
# Import Data
#......................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dtsrvy <- makecd2013survey(survey = dt)
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
#......................
# Subset to Pv
#......................
longlat <- dt %>% 
  dplyr::select(c("hv001", "longnum", "latnum")) %>% 
  dplyr::filter(!duplicated(.))

pvclust.weighted.nosf <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = srvyr::survey_total(count), 
                   plsmdn = srvyr::survey_total(pv18s, na.rm = T), 
                   plsmdprev = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95))

# need to keep integers, so will round
pvclust.weighted.nosf <- pvclust.weighted.nosf %>% 
  dplyr::mutate(plsmdn = round(plsmdn, 0),
                n = round(n, 0)) %>% 
  dplyr::left_join(., longlat, by = "hv001") # add back in longlat

#-------------------------------------------------------------------------
# Aggregate Covariates
#-------------------------------------------------------------------------
# precipitation already in dataset 
pvclst.covar <- dt %>% 
  dplyr::select("hv001", "precip_mean_cont_scale_clst") %>% 
  dplyr::filter(!duplicated(.))

# bring in cropland
crop <- readRDS("data/derived_data/vividepi_cropland_propmeans.rds") %>% 
  dplyr::select(c("hv001", "cropprop_cont_scale_clst"))
# bring in urban friction
urbanfric <- readRDS("data/derived_data/vividepi_fricurban_clstmeans.rds") %>% 
  dplyr::select(c("hv001", "frctmean_cont_scale_clst"))

pvclst.covar <- dplyr::left_join(x = pvclst.covar, y = crop, by = "hv001") %>% 
  dplyr::left_join(x = ., y = urbanfric, by = "hv001")

# join together
pvclust.weighted.nosf <- dplyr::left_join(pvclust.weighted.nosf, pvclst.covar, by = "hv001")

riskvars <- c("precip_mean_cont_scale_clst", 
              "cropprop_cont_scale_clst", "frctmean_cont_scale_clst")



#-------------------------------------------------------------------------
# Make Model Framework
#-------------------------------------------------------------------------
mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                formula = c("plsmdn ~ 1", 
                                            paste0("plsmdn ~ ", paste(riskvars, collapse = "+")))
)

coords <- as.formula(paste0("~ longnum + latnum"))
mod.framework$data <- lapply(1:nrow(mod.framework), function(x) return(pvclust.weighted.nosf))
mod.framework$coords <- lapply(1:nrow(mod.framework), function(x) return(coords))



####################################################################################
###########           PRIORS  & Direction Diagnostic                      ##########
####################################################################################
fit.glm <- glm(cbind(plsmdn, n - plsmdn) ~ 1, 
               data = pvclust.weighted.nosf,
               family = binomial)


mypriors.intercept <- PrevMap::control.prior(beta.mean = 0,
                                             beta.covar = 1,
                                             log.normal.nugget = c(0, 2.5), # this is tau2
                                             uniform.phi = c(0, 10),
                                             log.normal.sigma = c(0, 2.5))

# NB covar matrix
covarsmat <- matrix(0, ncol = 4, nrow=4) # three risk factors, 4 betas
diag(covarsmat) <- 1 # identity matrix

mypriors.mod <- PrevMap::control.prior(beta.mean = c(0, 0, 0, 0),
                                       beta.covar = covarsmat,
                                       log.normal.nugget = c(0, 2.5), # this is tau2
                                       uniform.phi = c(0, 10),
                                       log.normal.sigma = c(0, 2.5))

mcmcdirections.intercept <- PrevMap::control.mcmc.Bayes(burnin = 5e3, 
                                                        n.sim = 5e4+5e3,
                                                        thin = 10, 
                                                        L.S.lim = c(5,50),
                                                        epsilon.S.lim = c(0.01, 0.1),
                                                        start.nugget = 1,
                                                        start.sigma2 = 0.2,
                                                        start.beta = -5,
                                                        start.phi = 0.5,
                                                        start.S = predict(fit.glm))

mcmcdirections.mod <- PrevMap::control.mcmc.Bayes(burnin = 5e3, 
                                                  n.sim = 5e4+5e3,
                                                  thin = 10, 
                                                  L.S.lim = c(5,50),
                                                  epsilon.S.lim = c(0.01, 0.1),
                                                  start.nugget = 1,
                                                  start.sigma2 = 0.2, 
                                                  start.beta = c(-4, -0.2, 0, -0.02),
                                                  start.phi = 0.5,
                                                  start.S = predict(fit.glm))


mod.framework$mypriors <- lapply(1:nrow(mod.framework), 
                                 function(x) return(mypriors.intercept))
mod.framework$mcmcdirections <- lapply(1:nrow(mod.framework), 
                                       function(x) return(mcmcdirections.intercept))


# NOTE THIS HACK here, need to have multiple starting betas for the multiple betas in the model
mod.framework$mypriors[[ which(stringr::str_detect(mod.framework$formula, "\\+")) ]] <- mypriors.mod
mod.framework$mcmcdirections[[ which(stringr::str_detect(mod.framework$formula, "\\+")) ]] <- mcmcdirections.mod

# replicate this four times for our four chains
mod.framework <- lapply(1:4, function(x) return(mod.framework)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(name)


#......................
# Make a wrapper for PrevMap diagnostic chain
#......................
fit_bayesmap_wrapper <- function(path){
  input <- readRDS(path)
  name <- input$name  
  formula <- input$formula  
  trials = "n"
  coords <- input$coords[[1]]
  data <- input$data[[1]]
  mypriors <- input$mypriors[[1]]
  mcmcdirections <- input$mcmcdirections[[1]]
  kappa = 0.75
  
  ret <- PrevMap::binomial.logistic.Bayes(
    formula = as.formula(formula),
    units.m = as.formula(paste("~", trials)),
    coords = coords,
    data = data,
    control.prior = mypriors,
    control.mcmc = mcmcdirections,
    kappa = kappa
  )
  
  #......................
  # save out
  #......................
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_diagn_runs/", 
             recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_diagn_runs/",
                   name, ".diagnostic_run_ret.RDS")
  saveRDS(ret, file = outpath)
  return(0)
}




####################################################################################
###########            Diagnostic Run Drake Plan                           #########
#####################################################################################
setwd("analyses/07-spatial_prediction")
#......................
# make diagnostic map
#......................
# PrevMap is having trouble with the drake nesting. resorted to reading in input files
dir.create("prevmap_drake_inputs")
mod.framework$name <- paste0(mod.framework$name, seq(1, 4))
mod.framework.list <- split(mod.framework, 1:nrow(mod.framework))
names(mod.framework.list) <- purrr::map_chr(mod.framework.list, "name")
lapply(mod.framework.list, function(x){
  saveRDS(x, paste0("prevmap_drake_inputs/", x$name, "_prevmap_diag_input.RDS"))
})

diagruns <- list.files("prevmap_drake_inputs/",
                       full.names = TRUE)
diagruns <- tibble::tibble(path = diagruns)

# make drake map
diag_analysis_names <- stringr::str_split_fixed(basename(diagruns$path), "_", n = 2)[,1]
plan_diag <- drake::drake_plan(
  fits = target(
    fit_bayesmap_wrapper(path), 
    transform = map(
      .data = !!diagruns,
      .names = !!diag_analysis_names
    )
  )
)
####################################################################################
###########                          LONG Run                              #########
####################################################################################
# note that PrevMap::spatial.pred.binomial.Bayes access the "formula"
# call within the model object for prediction (with covariates). As a result,
# we can't use a wrapper with purrr as above. 
dir.create("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_long_runs/", 
           recursive = TRUE)
# Directions LONG RUN                      
mcmcdirections.intercept.long <- PrevMap::control.mcmc.Bayes(burnin = 5e4, 
                                                             n.sim = 5e5 + 5e4,
                                                             thin = 100, 
                                                             L.S.lim = c(5,50),
                                                             epsilon.S.lim = c(0.01, 0.1),
                                                             start.nugget = 1,
                                                             start.sigma2 = 0.2,
                                                             start.beta = -5,
                                                             start.phi = 0.5,
                                                             start.S = predict(fit.glm))

mcmcdirections.mod.long <- PrevMap::control.mcmc.Bayes(burnin = 5e4, 
                                                       n.sim = 5e5 + 5e4,
                                                       thin = 100, 
                                                       L.S.lim = c(5,50),
                                                       epsilon.S.lim = c(0.01, 0.1),
                                                       start.nugget = 1,
                                                       start.sigma2 = 0.2, 
                                                       start.beta = c(-4, -0.2, 0, -0.02),
                                                       start.phi = 0.5,
                                                       start.S = predict(fit.glm))



#......................
# drake plan for long intercept
#......................
plan_long_intercept <- drake::drake_plan(
  longrun_intercept = target(
    PrevMap::binomial.logistic.Bayes(
      formula = as.formula("plsmdn ~ 1"),
      units.m = as.formula("~ n"),
      coords = as.formula("~ longnum + latnum"),
      data = pvclust.weighted.nosf,
      control.prior = mypriors.intercept,
      control.mcmc = mcmcdirections.intercept.long,
      kappa = 0.75
    )),
  savelongrun_intercept = target(
    saveRDS(longrun_intercept, 
            file = "/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_long_runs/intercept_model.RDS"),
    hpc = FALSE
  )
)


#......................
# drake plan for long covar
#......................
plan_long_covarmad <- drake::drake_plan(
  longrun_covarmod = target(
    PrevMap::binomial.logistic.Bayes(
      formula = as.formula("plsmdn ~ 1 + precip_mean_cont_scale_clst + 
                       cropprop_cont_scale_clst + frctmean_cont_scale_clst"),
      units.m = as.formula("~ n"),
      coords = as.formula("~ longnum + latnum"),
      data = pvclust.weighted.nosf,
      control.prior = mypriors.mod,
      control.mcmc = mcmcdirections.mod.long,
      kappa = 0.75
    )),
  savelongrun_covar = target(
    saveRDS(longrun_covarmod, 
            file = "/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_long_runs/covariate_model.RDS"),
    hpc = FALSE))



#............................................................
# bring drake plans together
#...........................................................
plan <- drake::bind_plans(plan_diag, 
                          plan_long_intercept,
                          plan_long_covarmad)

#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_clst/slurm_clustermq_LL_fit.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = 10,
     log_make = "prevamp_fit_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")

