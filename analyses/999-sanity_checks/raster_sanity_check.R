#............................................................
# understanding rasters with +proj=longlat
# https://www.rdocumentation.org/packages/raster/versions/3.4-5/topics/extract
#...........................................................

library(raster)

r <- raster(ncol=36, nrow=18, vals=1:(18*36))

###############################
# extract values by cell number
###############################
extract(r, c(1:2, 10, 100))
s <- stack(r, sqrt(r), r/r)
extract(s, c(1, 10, 100), layer=2, n=2)

###############################
# extract values with points
###############################
xy <- cbind(-50, seq(-80, 80, by=20))
extract(r, xy)

sp <- SpatialPoints(xy)
extract(r, sp, method='bilinear')

# examples with a buffer
extract(r, xy[1:3,], buffer=1000000)
extract(r, xy[1:3,], buffer=1000000, fun=mean)

## illustrating the varying size of a buffer (expressed in meters) 
## on a longitude/latitude raster
crs(r)
z <- extract(r, xy, buffer=1000000)
s <- raster(r)
for (i in 1:length(z)) { s[z[[i]]] <- i }
plot(s)

## compare with raster that is not longitude/latitude
crs(r) <- "+proj=utm +zone=17" 
xy[,1] <- 50
z <- extract(r, xy, buffer=8)
for (i in 1:length(z)) { s[z[[i]]] <- i }
plot(s)
