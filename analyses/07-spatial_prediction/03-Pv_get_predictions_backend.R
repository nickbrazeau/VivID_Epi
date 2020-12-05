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
tol <- 1e-3
set.seed(48, "L'Ecuyer")

#......................
# Import Data
#......................
intercept <- readRDS("analyses/07-spatial_prediction/prevmap_long_runs/intercept_model.RDS")
covar <- readRDS("analyses/07-spatial_prediction/prevmap_long_runs/covariate_model.RDS")
# to do add covars if needed
gp.mod.framework <- tibble::tibble(name = c("intercept", "covars"),
                                   mcmc = list(intercept, covar))


#...............................
# sample coordinates for prediction surface 
#...............................
# boundaries for prediction
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred.coords <- splancs::gridpts(poly, xs=0.05, ys=0.05)
colnames(grid.pred.coords) <- c("longnum","latnum")
grid.pred.coords.df <- as.data.frame(grid.pred.coords)
# convert to sf tidy  
grid.pred.coords <- sf::st_as_sf(grid.pred.coords.df, coords = c("longnum", "latnum"),
                                 crs = "+init=epsg:4326")

#...............................
# extract covariate information for prediction surface 
#...............................
# Raster surfaces for risk factors
riskvars <- c("precip_mean_cont_scale_clst", 
              "cropprop_cont_scale_clst", "accmean_cont_scale_clst")

precipraster <- readRDS("data/derived_data/vividepi_precip_study_period_effsurface.rds") 
# transform like fitting
precipraster <- raster::scale(precipraster, center = TRUE, scale = TRUE)

# take mean here since we converted to binary
cropraster <- readRDS("data/derived_data/vividepi_cropland_surface.rds")
cropraster <- raster::aggregate(cropraster, fact = 18, fun = mean)
# transform like fitting
values(cropraster) <- my.scale(logit(values(cropraster), tol = tol))


# mean here under the assumption that acess is measured continuously across space
accraster <- readRDS("data/derived_data/vividepi_accessurban_surface.rds")
accraster <- raster::aggregate(accraster, fact = 6, fun = mean)
values(accraster) <- my.scale(log(values(accraster) + tol))

# stack 
predcovars <- raster::stack(precipraster, cropraster, accraster)
names(predcovars) <- riskvars

pred.df <- raster::extract(
  x = predcovars,
  y = sf::as_Spatial(grid.pred.coords),
  buffer = 6000,
  fun = mean,
  na.rm = T,
  sp = F)

pred.df <- as.data.frame(pred.df)
pred.df <- cbind.data.frame(grid.pred.coords.df, pred.df)

#...............................
# Bounds so we only do interpolation
# (i.e. coerce extrapolation values)
#...............................
# precip
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
precip <- dt %>% 
  dplyr::select("hv001", "precip_mean_cont_scale_clst") %>% 
  dplyr::filter(!duplicated(.))
max.precip <- max(precip$precip_mean_cont_scale_clst)
min.precip <- min(precip$precip_mean_cont_scale_clst)

pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst > max.precip] <- max.precip
pred.df$precip_mean_cont_scale_clst[pred.df$precip_mean_cont_scale_clst < min.precip] <- min.precip


# bring in cropland
crop <- readRDS("data/derived_data/vividepi_cropland_propmeans.rds") 
max.crop <- max(crop$cropprop_cont_scale_clst)
min.crop <- min(crop$cropprop_cont_scale_clst)
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst > max.crop] <- max.crop
pred.df$cropprop_cont_scale_clst[pred.df$cropprop_cont_scale_clst < min.crop] <- min.crop


# bring in acction surface
acc <- readRDS("data/derived_data/vividepi_accurban_clstmeans.rds") 
max.acc <- max(acc$accmean_cont_scale_clst)
min.acc <- min(acc$accmean_cont_scale_cls)
pred.df$accmean_cont_scale_clst[pred.df$accmean_cont_scale_clst > max.acc] <- max.acc
pred.df$accmean_cont_scale_clst[pred.df$accmean_cont_scale_clst < min.acc] <- min.acc


#..............................................................
# Convert to raster for sampling
#..............................................................
covar.rstr.pred <- raster::rasterFromXYZ(pred.df, 
                                         res = c(0.05, 0.05),
                                         crs = "+init=epsg:4326")

# save this out for later
saveRDS(covar.rstr.pred, file = "data/derived_data/vividepi_spatial_covar_feature_engineer.rds")

# down sample
covar.rstr.pred.downsmpl <- raster::sampleRandom(covar.rstr.pred, 
                                                 size = 2e4,
                                                 na.rm = T,
                                                 xy = T,
                                                 sp = F)

colnames(covar.rstr.pred.downsmpl)[1:2] <- c("longnum", "latnum")

#...............................
# Setup Map Dataframe  
#...............................
# set up grid.pred
gp.mod.framework$grid.pred <- list(covar.rstr.pred.downsmpl[,c("longnum", "latnum")],
                                   covar.rstr.pred.downsmpl[,c("longnum", "latnum")])

# set up predictors
gp.mod.framework$predictors <- list(NULL, covar.rstr.pred.downsmpl[,c("precip_mean_cont_scale_clst", 
                                                                      "cropprop_cont_scale_clst",
                                                                      "accmean_cont_scale_clst")])



# make wrapper
pred_PrevMap_bayes_wrapper <- function(mcmc, grid.pred, predictors){
  ret <- PrevMap::spatial.pred.binomial.Bayes(object = mcmc, 
                                              grid.pred = grid.pred,
                                              predictors = predictors, 
                                              type = "marginal", 
                                              scale.predictions = "prevalence",
                                              quantiles = NULL, 
                                              standard.error = TRUE, 
                                              thresholds = NULL)
  return(ret)
  
}

####################################################################################
###########                 Preds with Drake                               #########
#####################################################################################
paramsdf <- gp.mod.framework[, c("mcmc", "grid.pred", "predictors")]
dir.create("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_predictions/",
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
            file = "/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_predictions/intercept_predictions.RDS"),
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
            file = "/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi/analyses/07-spatial_prediction/prevmap_predictions/covar_predictions.RDS"),
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
        clustermq.template = "drake_clst/slurm_clustermq_LL_pred.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = 2,
     log_make = "prevamp_pred_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)






sjob <- rslurm::slurm_apply(f = pred_PrevMap_bayes_wrapper, 
                            params = paramsdf, 
                            jobname = 'Prevmap_predictions_large',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = "500g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "11-00:00:00"))




