# https://www.rdocumentation.org/packages/energy/versions/1.7-6/topics/distance%20correlation
# https://cran.r-project.org/web/packages/energy/energy.pdf
# https://projecteuclid.org/euclid.aoas/1267453933
# http://yunus.hacettepe.edu.tr/~iozkan/eco742/Brownian.html
# https://projecteuclid.org/download/pdfview_1/euclid.aos/1201012979
# https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
library(energy)

adj.dist <- function(x, e=2){
  # Double center distance matrix
  xd <- as.matrix(dist(x, p=e))
  ln <- dim(xd)[1]
  # keep colmeans and rowmeans as distance matrix may be asymmetric 
  # (based on some sim.)
  cm <- colMeans(xd)
  rm <- rowMeans(xd)
  # grand mean (for only now.)!!
  gm <- mean(xd)
  for (i in 1:ln){
    xd[i,] <- xd[i,] - rm[i]
    xd[,i] <- xd[,i] - cm[i]
  }
  return(xd + gm)
}    

distance_cor <- function(x,y,ind=2){
  x_adj <- adj.dist(x, e=ind)
  y_adj <- adj.dist(y, e=ind) 
  v_xy <- mean(x_adj  *  y_adj)
  v_xx <- mean(x_adj  *  x_adj)
  v_yy <- mean(y_adj  *  y_adj)
  # distance covariance
  d.cov <- sqrt(v_xy)
  d.cor <- 0
  if (sqrt(v_xx * v_yy) != 0) d.cor <- d.cov/(sqrt(sqrt(v_xx) * sqrt(v_yy)))
  tst <- dim(x_adj)[1]*v_xy 
  return(data.frame(dcv=d.cov, dcr=d.cor, d.varx= sqrt(v_xx), d.vary= sqrt(v_yy), test=tst))
}

start_time <- Sys.time()
distance_cor(m[1:1581,1], m[1:1581,2])
end_time <- Sys.time()

end_time - start_time

cat("now to energy package")
start_time <- Sys.time()
energy::dcor(m[1:1581,1], m[1:1581,2])
end_time <- Sys.time()

end_time - start_time

cat("done")

x <- iris[1:50, 1:4]
y <- iris[51:100, 1:4]
dcor(x, y)
distance_cor(x,y)
