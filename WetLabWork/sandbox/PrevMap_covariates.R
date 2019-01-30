
# PrevMap_covariates.R

# Author: Bob Verity
# Date: 2017-10-23

# Purpose:
# Demonstrates the PrevMap package, fitting a Gaussian Process model to transformed frequencies along with covariates.

# ------------------------------------------------------------------

# load bobFunctions package. If not already installed, this can be obtained from github via the devtools command install_github('bobverity/bobFunctions')
library(bobFunctions)

# other packages
library(PrevMap)
library(raster)
library(geoR)

# ------------------------------------------------------------------

set.seed(2)

# we want to run PrevMap with some covariates. Read in elevation data as an example covariate
topoData <- raster("~/Dropbox/Bob/Work/Collaborations/Meshnick DRC/Analyses/2017.08.15 Missing children/R scripts/PrevMap examples/raleigh-w.DEM")

# subset to smaller space
nCell <- 300
topoData  <- as.matrix(topoData)[1:nCell,1:nCell]
topoData <- topoData - median(topoData)

# get into long format
grid.pred <- data.frame(long=as.vector(row(topoData)/nCell), lat=as.vector(col(topoData)/nCell), elevation=as.vector(topoData))

# add another covariate that is just random noise. Hopefully the model will drop this out as non-significant
grid.pred$fake_covariate <- rnorm(nrow(grid.pred))

# define true prevalence at all points as a function of covariates plus some Perlin noise (i.e. smoothly varying noise)
topoData_noise <- 0.1*topoData + 5*perlinNoise(nCell,nCell,10,10)
y <- -2 + as.vector(topoData_noise)
grid.pred$prev_true <- 1/(1+exp(-y))

# set parameters of simulated data
k <- 100	# number of spatial clusters
n <- 20		# number of samples per cluster

# simulate observed data as subset of grid
dat <- as.data.frame(grid.pred[sample(nrow(grid.pred),k),])

# generate binomial observation data from underlying true prevalence. The true prevalence is assumed unknown in the real data - all we have access to are these binomial counts.
dat$count <- rbinom(k,n,prob=dat$prev_true)
dat$prev_obs <- dat$count/n

# logit transform observed frequencies. eps value ensures that transformed values are finite
eps <- 1e-3
dat$y <- log((dat$prev_obs+eps)/(1-dat$prev_obs+eps))

# plot simple map of elevation and cluster locations
win(2,2)
plot(raster(topoData), asp=NA, col=rev(bobMapCols(100)), main="relative elevation\n(real covariate)")
points(dat$long, dat$lat, cex=0.5)

# plot (unknown) true prevalence
plot(raster(matrix(grid.pred$prev_true,nCell)), asp=NA, zlim=c(0,1), col=rev(bobFireIce(100)), main="true prevalence")

# ------------------------------------------------------------------

# fit model by MLE without covariates
fit.MLE <- linear.model.MLE(formula=y~1, coords=~long+lat, data=dat, start.cov.pars=c(1,1), kappa=0.5)

summary(fit.MLE, log.cov.pars=FALSE)

# make spatial predictions
reps <- 1e2

pred.MLE <- spatial.pred.linear.MLE(fit.MLE, grid.pred=grid.pred[,c("long","lat")], scale.predictions="prevalence", n.sim.prev=reps, standard.errors=TRUE)

# extract main result in matrix
grid.pred.prev <- matrix(pred.MLE$prevalence$predictions,nCell)

# plot predicted prevalence surface
plot(raster(grid.pred.prev), asp=NA, zlim=c(0,1), col=rev(bobFireIce(100)), main="estimated prevalence\nwithout covariates")

# ------------------------------------------------------------------

# fit model by MLE with covariates
fit.MLE <- linear.model.MLE(formula=y~1+elevation+fake_covariate, coords=~long+lat, data=dat, start.cov.pars=c(1,1), kappa=0.5)

summary(fit.MLE, log.cov.pars=FALSE)

# make spatial predictions
reps <- 1e2

pred.MLE <- spatial.pred.linear.MLE(fit.MLE, grid.pred=grid.pred[,c("long","lat")], predictors=as.data.frame(grid.pred), scale.predictions="prevalence", n.sim.prev=reps, standard.errors=TRUE)

# extract main result in matrix
grid.pred.prev <- matrix(pred.MLE$prevalence$predictions,nCell)

# plot predicted prevalence surface
plot(raster(grid.pred.prev), asp=NA, zlim=c(0,1), col=rev(bobFireIce(100)), main="estimated prevalence\nwith covariates")
