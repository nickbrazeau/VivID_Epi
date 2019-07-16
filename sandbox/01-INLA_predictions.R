#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create a spatial prediction model
# using INLA 
# Using MAP's codefrom their publication (PMC )
# https://github.com/malaria-atlas-project/malariaAtlas
# as guidance as well as several INLA tutorials:
# http://www.r-inla.org/events/newtutorialsonspatialmodelsininla
# https://haakonbakka.bitbucket.io/btopic122.html
# https://ourcodingclub.github.io/2018/12/04/inla.html
#
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
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "internal_datamap_files/DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))
dtsrvy <- makecd2013survey(survey = dt)

space <- dt[,c("hv001", "longnum", "latnum")]
sf::st_geometry(space) <- NULL
space <- space[!duplicated(space), ]

#TODO extand past intercept model

#----------------------------------------------------------------------------------------------------
# Mixed Effects Model, Non-Spatial
#----------------------------------------------------------------------------------------------------
eq <- as.formula(   paste("pv18s ~ 1 + (1|hv001)")  ) 

mme.sat.nonsp <- lme4::glmer(eq,
                             family = binomial(link="logit"),
                             data = dt )

merTools::fastdisp(mme.sat.nonsp)
broom.mixed::tidy(mme.sat.nonsp, conf.int = T, exponentiate = T)









#----------------------------------------------------------------------------------------------------
# GEOSTATISTICAL MODELING with INLA
#----------------------------------------------------------------------------------------------------
# The basic steps for the spde spatial random effect are to estimate the Mattern Covariance: 
#   (1) Make a "mesh" that allows for faster Mattern Covariance calculation
#   (2) Make an "A-matrix" for actual observed data points that tracks w/ mesh
#   (3) Calculate Mattern Covariance from Mesh
#   (4) Project calculations onto observed data points = "A-Matrix" 


# (1) make mesh
mesh <- inla.mesh.2d(space[, c('longnum', 'latnum')], 
                     max.edge = c(0.8,2)) # TODO determine what these should be
                    
plot(mesh) # visualise the mesh

# (2) make A-matrix
A <- inla.spde.make.A(mesh = mesh, 
                      loc = as.matrix(space[, c('longnum', 'latnum')]))

# Calculate Mattern Covariance (e.g. make spde object)
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            alpha = 2, 
                            prior.range = c(2,0.01), 
                            prior.sigma = c(2.7,0.1)) #TODO priors?


# make INLA stack for inla model
# https://haakonbakka.bitbucket.io/btopic122.html
stack = inla.stack(tag='pv_est',
                   data=list(y=dt$pv18s), # response
                   
                   effects=list(
                     # - The Model Components
                     s=1:spde$n.spde, 
                     # - The first is 's' (for spatial)
                     data.frame(intercept=1, x=dt[,spatial_covar])
                     # - The second is all fixed effects
                     ),
                   A=list(A, 1)
                   # - First projector matrix is for 's'
                   # - second is for 'fixed effects'
                   # note, if you add on more random effects, need additional 1s here (see https://ourcodingclub.github.io/2018/12/04/inla.html)
)



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









