

smd <- function(dat, target){
  
  # find levels of binary treatment
  lvls <- levels( factor(unlist(dat[,target])) )
  
  # pull data for group 1
  dat.m1 <- dat[ unlist( dat[,target] ) == lvls[1] , ]
  
  # pull data for group 1
  dat.m2 <- dat[ unlist( dat[,target] ) == lvls[2] , ]
  
  mean.m1 <- mean(unlist( dat.m1[,"val"] ))
  mean.m2 <- mean(unlist( dat.m2[,"val"] ))
  sd.m1 <- sd(unlist( dat.m1[,"val"] ))
  sd.m2 <- sd(unlist( dat.m2[,"val"] ))
  
  smd <- (mean.m1 - mean.m2) / sqrt(sd.m1 + sd.m2)
  smd <- abs( smd )
  return(smd)
  
  
}

# First need to create a null distribution
make.null.distribution.binaryTx <- function(target, covars, data){
  n <- nrow(data)
  # TODO reflected normal
  wi <- abs( rnorm(n, mean = 1, sd = 0.25) )
  
  
  
  data.map <- data %>% 
    dplyr::select(c(target, covars)) %>% 
    tidyr::gather(., key = "covar", value = "val", covars) %>% 
    dplyr::group_by(covar) %>% 
    tidyr::nest()
  
  data.map$widata <- purrr::map(data.map$data, function(x){
    x <- x[sample(1:nrow(x), size = n, prob = wi, replace = T), ]
    return(x)
  })
  
  data.map$smd <- purrr::map(data.map$widata, smd, target = target)
  
  avgdj <- mean( unlist(data.map$smd) )
  
  return(avgdj)
  
}




my.pearson <- function(matrix){
  ret <- cor(matrix, method = "pearson")
  ret <- ret[2,1] # give off diagonal
  return(ret)
}


# First need to create a null distribution
make.null.distribution.continuousTx <- function(target, covars, data){
  n <- nrow(data)
  # TODO reflected normal
  wi <- abs( rnorm(n, mean = 1, sd = 0.25) )
  
  
  data.map <- data %>% 
    dplyr::select(c(target, covars)) %>% 
    tidyr::gather(., key = "covar", value = "val", covars) %>% 
    dplyr::group_by(covar) %>% 
    tidyr::nest()
  
  data.map$widata <- purrr::map(data.map$data, function(x){
    x <- x[sample(1:nrow(x), size = n, prob = wi, replace = T), ]
    return(x)
  })
  
  data.map$cor <- purrr::map(data.map$widata, my.pearson)
  
  avgdj <- mean( unlist(data.map$cor) )
  
  return(avgdj)
  
  
}
