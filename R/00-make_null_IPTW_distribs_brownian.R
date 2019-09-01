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
  wi <- abs( rnorm(n, mean = 1, sd = 0.25) )
  # treating the weights here as your probability of being sampled from the gen pop
  # this is true for the long-run average (eg over many iterations)
  wirows <- sample(1:nrow(x), size = n, prob = wi, replace = T)
  data.wirows <- data[wirows, ]
  
  data.list <- lapply(1:length(covars), 
                      function(x){return(as.data.frame(data.wirows[, c(target, covars[x]) ]))})
  

  data.dist <- lapply(data.list, function(x){
    
    if(is.factor(x[,1])){ # note, only have binary factors so this ok
      x[,1] <- as.numeric(x[,1])
    }
    
    if(is.factor(x[,2])){
      x[,2] <- as.numeric(x[,2])
    }
    
    ret <- energy::dcor(x = x[, 1], y = x[, 2])
    
    return(ret)
  })
 
  # make final return
  dij <- unlist(data.dist)
  dij <- mean(dij)
  return(dij)
  
}
  
  