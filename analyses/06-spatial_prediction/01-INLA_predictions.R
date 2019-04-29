#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create a spatial prediction model
# using INLA 
# Using MAP's codefrom their publication (PMC )
# https://github.com/malaria-atlas-project/malariaAtlas
# as guidance
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(srvyr) #wrap the survey package in dplyr syntax
library(INLA)

#......................
# Import Data
#......................
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")
pvprev <- mp %>% 
  dplyr::filter(maplvl == "hv001" & plsmdmspec == "pv18s") %>% 
  tidyr::unnest(.)


#----------------------------------------------------------------------------------------------------
# GEOSTATISTICAL MODELING with INLA
#----------------------------------------------------------------------------------------------------

# Building a mesh & spatial weights matrix

mesh <- inla.mesh.2d(pvprev[, c('longnum', 'latnum')], 
                     max.edge = c(0.8,2), 
                     cutoff = 0, 
                     min.angle = 21, 
                     offset = c(2, 5)) 
plot(mesh) # visualise the mesh

# Build weight matrix
A <- inla.spde.make.A(mesh = mesh, 
                      loc = as.matrix(pvprev[, c('longnum', 'latnum')]))

# Create spde2 object
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, 
                            prior.range = c(2,0.01), 
                            prior.sigma = c(2.7,0.1))

# stack spatial objects and data.
stk.acc <- inla.stack(tag = 'estimation', 
                      data = list(pos = pvpr$pv_pos_1to99),
                      A = list(A, 1),
                      effects = list(space = 1:spde$n.spde, 
                                     data.frame(B0 = 1, 
                                                pvpr[,c("NightTimeTemp.scaled", 
                                                        "logElevation.scaled",
                                                        "Rainfall.scaled",
                                                        "logTravelToCities.scaled")]))
                      )

# define model formula and run INLA model
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









