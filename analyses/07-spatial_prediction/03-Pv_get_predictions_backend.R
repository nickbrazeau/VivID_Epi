#----------------------------------------------------------------------------------------------------
# Purpose of this script is to get spatial predictions from 
# spatial models produced from prevmap
#----------------------------------------------------------------------------------------------------
source("R/00-functions_maps.R") 
source("R/00-functions_basic.R")
library(tidyverse)
library(drake)
library(raster)
library(PrevMap)
set.seed(48, "L'Ecuyer")


#............................................................
# set up spatial model
#...........................................................
#......................
# Import Data
#......................
intercept <- readRDS("analyses/07-spatial_prediction/prevmap_long_runs/intercept_model.RDS")
covar <- readRDS("analyses/07-spatial_prediction/prevmap_long_runs/covariate_model.RDS")
# to do add covars if needed
gp.mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                   mcmc = list(intercept, covar))


#......................
# read in prediction df and downsample for spatial predictions
#......................
# read in spatial raster preds
covar.rstr.pred <- readRDS("data/derived_data/vividepi_spatial_covar_feature_engineer.rds")

# convert from raster to array/matrix
covar.rstr.pred.downsmpl <- raster::sampleRandom(covar.rstr.pred,
                                                 size = raster::ncell(covar.rstr.pred), # want  all locations 
                                                 na.rm = T,
                                                 xy = T,
                                                 sp = F
                                                 )
# sanity check
sum(duplicated(covar.rstr.pred.downsmpl))


# prevamp needs this as a dataframe
# also need to rename x and y names
colnames(covar.rstr.pred.downsmpl)[1:2] <- c("longnum", "latnum")
covar.rstr.pred.downsmpl <- as.data.frame(covar.rstr.pred.downsmpl)

#...............................
# bring together Map Dataframe  
#...............................
# bring in grid.pred
gp.mod.framework$grid.pred <- list(covar.rstr.pred.downsmpl[,c("longnum", "latnum")],
                                   covar.rstr.pred.downsmpl[,c("longnum", "latnum")])

# bring in predictors
gp.mod.framework$predictors <- list(NULL, covar.rstr.pred.downsmpl[,c("precip_mean_cont_scale_clst", 
                                                                      "hlthdist_cont_scale_clst")])




####################################################################################
###########                 Preds with Drake                               #########
#####################################################################################
paramsdf <- gp.mod.framework[, c("mcmc", "grid.pred", "predictors")]
dir.create("analyses/07-spatial_prediction/prevmap_predictions/",
           recursive = TRUE)

#......................
# drake plan for prediction of intercept
#......................
plan_pred_intercept <- drake::drake_plan(
  pred_intercept = target(
    PrevMap::spatial.pred.binomial.Bayes(object = paramsdf$mcmc[[1]], 
                                         grid.pred = paramsdf$grid.pred[[1]],
                                         predictors = paramsdf$predictors[[1]], 
                                         type = "marginal", 
                                         scale.predictions = "prevalence",
                                         quantiles = NULL, 
                                         standard.error = TRUE, 
                                         thresholds = NULL)
    ),
  save_pred_intercept = target(
    saveRDS(pred_intercept, 
            file = "analyses/07-spatial_prediction/prevmap_predictions/intercept_predictions.RDS"),
    hpc = FALSE
  )
)


#......................
# drake plan for prediction of covariate
#......................
plan_pred_covar <- drake::drake_plan(
  pred_covar = target(
    PrevMap::spatial.pred.binomial.Bayes(object = paramsdf$mcmc[[2]], 
                                         grid.pred = paramsdf$grid.pred[[2]],
                                         predictors = paramsdf$predictors[[2]], 
                                         type = "marginal", 
                                         scale.predictions = "prevalence",
                                         quantiles = NULL, 
                                         standard.error = TRUE, 
                                         thresholds = NULL)
  ),
  save_pred_covar = target(
    saveRDS(pred_covar, 
            file = "analyses/07-spatial_prediction/prevmap_predictions/covar_predictions.RDS"),
    hpc = FALSE
  )
)

#......................
# bring together
#......................
plan <- drake::bind_plans(plan_pred_intercept,
                          plan_pred_covar)

#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "analyses/07-spatial_prediction/drake_clst/slurm_clustermq_LL_pred.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = 2,
     log_make = "analyses/07-spatial_prediction/prevmap_pred_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, 
     lock_cache = FALSE)






