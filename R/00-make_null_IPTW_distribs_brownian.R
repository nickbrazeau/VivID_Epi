#--------------------------------------------------------------------------
# Purpose of this script is to make a null distribution
# https://www.rdocumentation.org/packages/energy/versions/1.7-6/topics/distance%20correlation
# https://cran.r-project.org/web/packages/energy/energy.pdf
# https://projecteuclid.org/euclid.aoas/1267453933
# http://yunus.hacettepe.edu.tr/~iozkan/eco742/Brownian.html
# https://projecteuclid.org/download/pdfview_1/euclid.aos/1201012979
# https://blogs.sas.com/content/iml/2018/04/04/distance-correlation.html
#--------------------------------------------------------------------------


make.null.distribution.energy <- function(target, data, covars){
  n <- nrow(data)
  # TODO reflected normal
  wi <- abs( rnorm(n, mean = 1, sd = 0.25) )
  
  
  data.list <- lapply(1:length(covars), 
                      function(x){return(as.data.frame(data[, c(target, covars[x]) ]))})
  
  data.list.rand <- lapply(data.list, function(x){
    x <- x[sample(1:nrow(x), size = n, prob = wi, replace = T), ]
    return(x)
  })
  
  data.dist <- lapply(data.list.rand, function(x){
    
    if(is.factor(x[,1])){
      x[,1] <- as.numeric(x[,1])
    }
    
    if(is.factor(x[,2])){
      x[,2] <- as.numeric(x[,2])
    }
    
    nrand <- floor(nrow(x) * 0.1) # take a 10% sample for speed
    nrand.rows <- sample(x = 1:nrow(x), size = nrand, replace = F) # pull random rows
    ret <- energy::dcor(x = x[nrand.rows, 1], y = x[nrand.rows, 2])
    return(ret)
  })
 
  # make final return
  dij <- unlist(data.dist)
  dij <- mean(dij)
  return(dij)
  
   
}
  
  