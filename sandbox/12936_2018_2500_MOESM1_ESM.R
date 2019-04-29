# from their publication PMC
# https://github.com/malaria-atlas-project/malariaAtlas
### Additional file 1: Mock Analysis code ####
#### Mock Analysis 1: ####
## Predicting the spatial distribution of Plasmodium vivax using malariaAtlas-derived response and covariate data

## Loading packages
# Data access, management and visualisation
library(malariaAtlas)
library(boot)
library(ggplot2)

# Spatial utility and modelling
library(INLA)
library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(seegSDM) # available via devtools::install_github("SEEG-Oxford/seegSDM")
library(seegMBG) # available via devtools::install_github("SEEG-Oxford/seegMBG")

## Download & Visualise Response data 

# Define spatial extent and download points within this area
extent <- matrix(c(-65.2,-11.8, -52, 1), 2, 2, dimnames = list(c("x", "y"), c("min", "max")))
pvpr_raw <- getPR(extent = extent, species = "pv")
pvpr_raw <- pvpr_raw[!is.na(pvpr_raw$pr),] # Subset to data points for which all data is publicly available  

# Visualise these points 
shp <-  as(raster::extent(extent), "SpatialPolygons") # create a shapefile of arbitrary study extent
shp_df <- ggplot2::fortify(shp) # convert to a data.frame for ggplot-mapping

p <- autoplot(pvpr_raw,
              shp_df = shp_df,
              printed = FALSE,
              map_title = "PvPR Surveys in Example Study Area", 
              fill_legend_title = "Raw PvPR")
print(p)

# Use convertPrevalence to standarize pvpr to Pv parasite rate in age range 1 - 99 
pvpr_raw$pv_pr_1to99 <- convertPrevalence(prevalence = pvpr_raw$pr, 
                                          age_min_in = pvpr_raw$lower_age, 
                                          age_max_in =  pvpr_raw$upper_age,
                                          age_min_out = rep(1, length(pvpr_raw$pr)),
                                          age_max_out = rep(99, length(pvpr_raw$pr)),
                                          parameters = "Pv_Gething2012")

# Replace points layer in plot 'p' with age-standardized points & visualise
p$layers[[2]] <-  geom_point(data = pvpr_raw, aes(x = longitude, y = latitude, fill = pv_pr_1to99, size = examined), shape = 21) 
p <- p + scale_fill_distiller(name = "PvPR1-99", palette = "RdYlBu")

print(p)

## Prepare Covariates  

# load in environmental covariates pre-cropped to study region
covariates <- raster::stack(raster::raster("LST_night.tif"), # NASA LP DAAC. Land surface temperature and emissivity 8-day L3 global 1km. version 005.  2015. https://lpdaac.usgs.gov. Accessed Feb 2017.
                            raster::raster("elevation.tif"), # NASA LP DAAC. SRTMGL3S: NASA Shuttle Radar Topography Mission Global 3 arc second sub-sampled. Version 003.  2013. https://lpdaac.usgs.gov. Accessed Mar 2016.
                            raster::raster("rainfall.tif")) # Hijmans RJ, Cameron SE, Parra JL, Jones PG, Jarvis A. Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology. 2005; 25:1965-1978.

# download raster of travel time to cities (Weiss et al 2018) for study area & visualise this
access <- malariaAtlas::getRaster(surface = "A global map of travel time to cities to assess inequalities in accessibility in 2015",
                                  extent = extent)
autoplot_MAPraster(access)

# stack all covariate rasters, rename these and visualise
covariates <- raster::stack(covariates, access)
names(covariates) <- c("NightTimeTemp", "Elevation", "Rainfall", "TravelToCities")
covariate_plot <- autoplot_MAPraster(covariates)

# Convert age-standardised PR value to age-standardised number of positive diagnoses to enable weighting by sample size while using binomial likelihood in inla model
pvpr_raw$pv_pos_1to99 <- pvpr_raw$pv_pr_1to99*pvpr_raw$examined

# apply transformations to covariates and cap extreme values (far beyond those in regions where we have response data)
covariates_trans <- covariates

# log-transform Elevation and Travel time to cities
covariates_trans$Elevation <- log(covariates$Elevation+22.001)
covariates_trans$TravelToCities <- log(covariates$TravelToCities+0.1)

# cap elevation and travel time to cities to range where we have response data
elev_cap <- c(min(pvpr_cov$Elevation), max(pvpr_cov$Elevation))
values(covariates_trans$Elevation)[values(covariates_trans$Elevation)>elev_cap[2]] <- elev_cap[2]
values(covariates_trans$Elevation)[values(covariates_trans$Elevation)<elev_cap[1]] <- elev_cap[1]

access_cap <- c(min(pvpr_cov$TravelToCities), max(pvpr_cov$TravelToCities))
values(covariates_trans$TravelToCities)[values(covariates_trans$TravelToCities)>access_cap[2]] <- access_cap[2]
values(covariates_trans$TravelToCities)[values(covariates_trans$TravelToCities)<access_cap[1]] <- access_cap[1]

# Scale covariates
covariates_trans <- raster::scale(covariates_trans)
names(covariates_trans) <- c("NightTimeTemp.scaled", "logElevation.scaled", "Rainfall.scaled", "logTravelToCities.scaled")

# extract covariate values for data point locations
pvpr <- cbind(pvpr_raw, raster::extract(covariates_trans, as.matrix(pvpr_raw[,c("longitude", "latitude")])))

## Geostatistical Model with INLA 

# Building a mesh & spatial weights matrix

mesh <- inla.mesh.2d(pvpr[, c('longitude', 'latitude')], 
                     max.edge = c(0.8,2), 
                     cutoff = 0, 
                     min.angle = 21, 
                     offset = c(2, 5)) 
plot(mesh) # visualise the mesh

# Build weight matrix
A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pvpr[, c('longitude', 'latitude')]))

# Create spde2 object
# Preliminary analysis of the variance of the residuals from non-spatial 
# 'base' models suggested a reasonable penalised complexity prior for sigma of
#  prior.sigma = c(7, 0.1)
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, prior.range = c(2,0.01), prior.sigma = c(2.7,0.1))

# stack spatial objects and data.
stk.acc <- inla.stack(tag = 'estimation', 
                      data = list(pos = pvpr$pv_pos_1to99),
                      A = list(A, 1),
                      effects = list(space = 1:spde$n.spde, 
                                     data.frame(B0 = 1, pvpr[,c("NightTimeTemp.scaled", "logElevation.scaled","Rainfall.scaled","logTravelToCities.scaled")]))
)

# %>% %>% %>% %>% %>% define model formula and run INLA model
formula.acc <- pos ~ -1 + B0 +  NightTimeTemp.scaled + logElevation.scaled + Rainfall.scaled+ logTravelToCities.scaled+ f(space, model = spde)## RF term
model.acc <- inla(formula.acc, 
                  family = 'binomial',
                  Ntrials = pvpr$examined,
                  data = inla.stack.data(stk.acc), 
                  control.predictor = list(A = inla.stack.A(stk.acc)), 
                  control.compute = list(waic = TRUE))

# look at model summaries
plot(model.acc, CI = TRUE)
summary(model.acc)

# predict pvpr over study area using Model 2
pred.acc <- predictRasterINLA(model.acc, raster = covariates_trans, 
                              mesh = mesh, type = "response", 
                              method = "MAP", constants = list("B0" = 1))

names(pred.acc) <- "INLA prediction w/ access"

# visualise predictions using malariaAtlas autoplot methods and ggplot2 adjustments
plot_acc_log <- autoplot_MAPraster(pred.acc, printed = FALSE)[[1]] +
  scale_fill_distiller(name = "PvPR1-99\n",palette = "RdYlBu", trans = "log10", limits = c(0.001, NA), breaks = c(0.001, 0.01, 0.1, 0.5), labels = c("Low", "", "", "High"), oob = scales::squish)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ggtitle("Predicted PvPR")

#### Mock Analysis 2: ####
## Testing a new modelling approach using in-built malariaAtlas zoon modules 

## Loading Packages 
library(zoon)
library(dplyr)

set.seed(26022018)

## Define and run zoon workflow for models WITHOUT mosquito occurrence data 

w <- workflow(Chain(malariaAtlas_PR(extent = c(44, 49, -20, -15), ISO = NULL, species = 'Pf', fold = 0),# Getting response data
                    malariaAtlas_PR(extent = c(44, 49, -24, -20), ISO = NULL, species = 'Pf', fold = 1)),
              Chain(Bioclim(extent = c(44, 49, -24, -15), resolution = 2.5, layers = c(1, 4, 12, 15)), # Getting covariate data
                    malariaAtlas_covariates(surface = 'A global map of travel time to cities to assess inequalities in accessibility in 2015',
                                            extent = c(44, 49, -24, -15))),
              RemoveNAs,
              LogisticRegression,
              Chain(PrintOccurrenceMap, 
                    PrintMap, 
                    PerformanceMeasures))

## Define and run zoon workflow for models WITH mosquito occurrence data 

w2 <- workflow(Chain(malariaAtlas_PR(extent = c(44, 49, -20, -15), ISO = NULL, species = 'Pf', fold = 0),# Getting response data
                     malariaAtlas_PR(extent = c(44, 49, -24, -20), ISO = NULL, species = 'Pf', fold = 1),
                     SpOcc(species = 'Anopheles arabiensis', 
                           extent =   c(44, 49, -24, -20)),
                     SpOcc(species = 'Anopheles gambiae', 
                           extent =   c(44, 49, -24, -20))),
               Chain(Bioclim(extent = c(44, 49, -24, -15), resolution = 2.5, layers = c(1, 4, 12, 15)), # Getting covariate data
                     malariaAtlas_covariates(surface = 'A global map of travel time to cities to assess inequalities in accessibility in 2015',
                                             extent = c(44, 49, -24, -15))),
               RemoveNAs,
               LogisticRegression,
               Chain(PrintOccurrenceMap, 
                     PrintMap, 
                     PerformanceMeasures))

## Summarise results and compare AUC values 

# count individuals by infection status
w$model.output[[1]]$data %>% group_by(fold, value) %>% count()
w2$model.output[[1]]$data %>% group_by(fold, value) %>% count()

# count individuals
w$model.output[[1]]$data %>% group_by(fold) %>% count()
w2$model.output[[1]]$data %>% group_by(fold) %>% count()

# count locations
w$model.output[[1]]$data %>% distinct(latitude, longitude, .keep_all = TRUE) %>% group_by(fold) %>% count()
w2$model.output[[1]]$data %>% distinct(latitude, longitude, .keep_all = TRUE) %>% group_by(fold) %>% count()

# retrieve AUC values from the model(s)
w$report[[1]][[3]]$auc
w2$report[[1]][[3]]$auc

## Convert raster outputs to MAPraster format for visualisation
w_df <- as.MAPraster(w$report[[1]][[2]])
w2_df <- as.MAPraster(w2$report[[1]][[2]])

w_plot <- autoplot(w_df, plot_titles = "", legend_title = "PfPR", printed = FALSE)[[1]] +
  scale_fill_distiller(name = "PfPR2-10\n",palette = "Blues", limits = c(0, 1), direction = 1,labels = c("Low", "","","","High"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 15))+
  theme(legend.text = element_text(size = 16), legend.title = element_text(face = "bold", size = 17), strip.text = element_text(face = "bold", size = 20))

w2_plot <- autoplot(w2_df, plot_titles = "", legend_title = "PfPR", printed = FALSE)[[1]] +
  scale_fill_distiller(name = "PfPR2-10\n",palette = "Blues", limits = c(0, 1), direction = 1, labels = c("Low", "","","","High"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 15))+
  theme(legend.text = element_text(size = 16), legend.title = element_text(face = "bold", size = 17), strip.text = element_text(face = "bold", size = 20))

