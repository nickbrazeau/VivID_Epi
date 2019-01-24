#--------------------------------------------------------
# Goal of this script is to simulate a process that depends 
# on space and has baseline confounding
#--------------------------------------------------------
# https://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise
# https://livefreeordichotomize.com/2019/01/17/understanding-propensity-score-weighting
library(tidyverse)
library(raster)
library(sp)
source("~/Documents/GitHub/VivID_Epi/analyses/00-functions_guassmap.R") 
set.seed(44)

#--------------------------------------------------------
# Spatial Dependencies
#--------------------------------------------------------
perlin_noise <- function( 
  n = 5,   m = 7,    # Size of the grid for the vector field
  N = 100, M = 100   # Dimension of the image
) {
  # For each point on this n*m grid, choose a unit 1 vector
  vector_field <- apply(
    array( rnorm( 2 * n * m ), dim = c(2,n,m) ),
    2:3,
    function(u) u / sqrt(sum(u^2))
  )
  f <- function(x,y) {
    # Find the grid cell in which the point (x,y) is
    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < n || j < m )
    # The 4 vectors, from the vector field, at the vertices of the square
    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]
    # Vectors from the point to the vertices
    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
    # Scalar products
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )
    # Weighted average of the scalar products
    s <- function(p) 3 * p^2 - 2 * p^3
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  xs <- seq(from = 1, to = n, length = N+1)[-(N+1)]
  ys <- seq(from = 1, to = m, length = M+1)[-(M+1)]
  outer( xs, ys, Vectorize(f) )
}

# create space probabilities based on perline noise (smoothly varying)
PN <- perlin_noise(n=2, m=3, N=10, M=10)
image(PN)

plot(raster(PN))

s <- matrix(0.1, 50, 100)
t <- matrix(1, 50, 100)
PN <- rbind(s,t)
plot(raster(PN))
#--------------------------------------------------------
# Setup Covariates
#--------------------------------------------------------
n <- 1e3

X <- mvtnorm::rmvnorm(n,
                      mean = c(0.5, 1),
                      sigma = matrix(c(2, 1, 1, 1), ncol = 2))


# set parameters of simulated data
k <- 100	# number of spatial points
m <- n/k		# number of samples per points

# simulate observed data as subset of grid
# get into long format
grid.pred <- data.frame(long=as.vector(row(PN)/nrow(PN)), 
                        lat=as.vector(col(PN)/nrow(PN)), pn=as.vector(PN))

# simulate observed data as subset of grid
clst <- as.data.frame(grid.pred[sample(nrow(grid.pred),k),])
clst$hv001 <- seq(1:nrow(clst))
clst <- lapply(1:(m/2), function(x){rbind(clst, clst)}) %>% 
  dplyr::bind_rows()

#--------------------------------------------------------
# SIMDATA
#--------------------------------------------------------
dat <- tibble(
  x_1 = X[, 1],
  x_2 = X[, 2],
  lon = unlist(clst[,1]),
  lat = unlist(clst[,2]),
  pn = unlist(clst[,3]),
  hv001 = unlist(clst[,4]),
  treatment = as.numeric(- 0.5 + 0.25 * x_1 + 0.75 * x_2 + rnorm(n, 0, 1) + pn*5 > 0),
  outcome = as.numeric(rbernoulli(n, plogis(2 * treatment + rnorm(n, 0, 1))))
)

PN ISN'T GET PICKED UP BC COVAR IS TOO MUCH FROM X1 AND X2'


#--------------------------------------------------------
# Check output
#--------------------------------------------------------

m1 <- glm(formula = outcome ~ treatment,
    data = dat,
    family = binomial(link = "logit"))
broom::tidy(m1, exponentiate = T, conf.int = T)
m2 <- glm(formula = outcome ~ treatment + x_1 + x_2,
    data = dat,
    family = binomial(link = "logit"))
broom::tidy(m2, exponentiate = T, conf.int = T)


# make cluster level
mp <- dat %>% 
  group_by(hv001, lat, lon) %>% 
  summarise(n=n(),
            outprev = mean(outcome))


# prevmap out and get a largely predict grid via long format
grid.pred <- expand.grid(list(lon=seq(0,1, by = 0.01), 
                        lat=seq(0,1, by = 0.01)))
prevraster <- fit_pred_spMLE(data = mp, outcome = "outprev", covar = "1", 
                      long_var = "lon", lat_var = "lat",
                      grid.pred = grid.pred[, c("lon", "lat")], kappa = 0.5, 
                      start.cov.pars = c(1,1),
                      pred.reps = 1)


grid.pred$prev <- prevraster$pred$prevalence$predictions
plot( raster::rasterFromXYZ(cbind(grid.pred[,1],
                            grid.pred[,2],
                            grid.pred[,3])) )
ggplot() + 
  geom_raster(data = grid.pred, aes(lon, lat, fill = prev), alpha = 0.8) +
  scale_fill_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = 0.67) 


ggplot() + 
  geom_point(data = mp, aes(lon, lat, color = outprev), alpha = 0.8) +
  scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = 0.67) 
plot(raster(PN))

