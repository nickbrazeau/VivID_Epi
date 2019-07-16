#--------------------------------------------------------------------------
# Purpose of this script is to make a null distribution
# https://www.rdocumentation.org/packages/energy/versions/1.7-6/topics/distance%20correlation
# https://cran.r-project.org/web/packages/energy/energy.pdf
# https://projecteuclid.org/euclid.aoas/1267453933
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
    
    ret <- energy::dcor(x = x[,1], y = x[,2])
    return(ret)
  })
 
  # make final return
  dij <- unlist(data.dist)
  dij <- mean(dij)
  return(dij)
  
   
}
  
  