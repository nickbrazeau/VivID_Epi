#--------------------------------------------------------
# Goal of this script is to simulate a process that depends 
# on space and has baseline confounding
#--------------------------------------------------------
library(tidyverse)
library(raster)
library(sp)

# TODO check out randomfields instead 

#--------------------------------------------------------
# Spatial Dependencies
#--------------------------------------------------------
perlin_noise <- function( 
  # function copied from
  # https://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise
  
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


#'@param mean.n.clust numeric; mean number of individuals in a cluster \(drawn from a gaussian distribution)
#'@param sd.n.clust numeric; sd for number of individuals \(drawn from a gaussian distribution)
#'@param nclusts numeric; number of clusters to draw from the \class{Raster} Matrix
#'@param PN.matrix matrix; a matrix of probabilities to that predict an underlying spatial process. Originally drawn from a Perlin Noise process
#'@param buffer numeric; the radius around the centroid cluster point that individual probabilites are being drawn. This is a simplified version of the DHS aggregation strategy

extract.PN_grid_values <- function(mean.n.clust = 50, sd.n.clust = 10, nclusts = 10, PN.matrix = NULL, buffer = 0.1){
  
  pn.raster <- raster::raster(PN.matrix)
  projection(pn.raster) <- "+proj=utm" # make this have long-lat so buffer behaves in extract
  longmin <- extent(pn.raster)@xmin
  longmax <- extent(pn.raster)@xmax
  
  latmin <- extent(pn.raster)@ymin
  latmax <- extent(pn.raster)@ymax
  
  # find cluster locations ("long/lat")
  coords <- data.frame(long = runif(n = nclusts, longmin, longmax),
                       lat = runif(n = nclusts, latmin, latmax)
                       )
  # draw probability masses from PN grid corresponding to cluster coordinates and buffer zone
  pn.points <- raster::extract(x = pn.raster, y = coords, buffer = buffer) # this is a buffer zone that we are pulling points around
  coords <- data.frame(cluster = as.character( 1:length(pn.points) ), coords, stringsAsFactors = F)
  names(pn.points) <- as.character( 1:length(pn.points) ) # for dplyr 
    
  # find the number of individuals in each cluster
  
  n.smpls <- ceiling( rnorm(n = nclusts, mean = mean.n.clust, sd = sd.n.clust) )
  
  # TODO error handle if n.smpls exceeds the number of pn.points 

  pn.points.df <- lapply(1:length(pn.points), function(x){
    ret <- data.frame(cluster = as.character( x ),
                      sprob = unlist(pn.points[[x]]),
                      stringsAsFactors = F)
    ret <- ret[ sample(1:nrow(ret), size = n.smpls[x], replace = F), ]
    
    return(ret)
  }) %>% 
    dplyr::bind_rows()
  
  
  out <- dplyr::left_join(coords, pn.points.df, by = "cluster")
  return(out)

  
}

# make spatial data
# create space probabilities based on perline noise (smoothly varying)
PN <- perlin_noise(n=2, m=3, N=100, M=100)

plot(raster(PN)) # going to treat the output of this as a probability mass
data <- extract.PN_grid_values(mean.n.clust = 50, sd.n.clust = 10,
                       PN.matrix = PN, buffer = 0.1)

# now make our covariates dependent on cluster and take space into account for whether or not you get outcome
data.list <- split(data, factor(data$cluster))
data.list <- lapply(data.list, function(dat){
  n = nrow(dat)
  ZX <- mvtnorm::rmvnorm(n, 
                         mean = c(0.5, 1),
                         sigma = matrix(c(2, 1, 1, 1), ncol = 2))
  dat <- tibble(
    z_1 = ZX[,1],
    z_2 = ZX[,2]
  ) %>% 
    dplyr::bind_cols(dat) %>% 
    dplyr::mutate(
      treatment = as.numeric(- 0.5 + 0.25 * z_1 + 0.33 * z_2  + rnorm(n, 0, 1) > 0),
      outcome = 2 * treatment + 1.5 *sprob + rnorm(n, 0, 1)
    )
  return(dat)

}) 


SimData <- data.list %>% dplyr::bind_rows(.)


plot(raster(PN)) # going to treat the output of this as a probability mas
points(SimData$lat, SimData$long)





