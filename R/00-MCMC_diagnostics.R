
#---------------------------------------------------
# Internal use function, not good for corner cases
#----------------------------------------------------

make_mcmc_chain_plots.carbayes.diag <- function(chaindat, filename){
    
    jpeg(filename, height = 8, width = 11, units = "in", res=200)
    par(mfrow=c(2,2))
    plot(chaindat[[1]][[1]])
    plot(chaindat[[1]][[2]])
    plot(chaindat[[1]][[3]])
    plot(chaindat[[1]][[4]])
    graphics.off()
  
}

wrap_chain_plotter.carbayes.diag <- function(diag.dir, chains){
  # this function does not return anything
  # it is internally making plots
  
  purrr::pmap(chains, function(data, name){
    
    filename <- paste0(diag.dir, name, "_", colnames(data), ".jpg") 
    
    for(i in 1:ncol(data)){
      make_mcmc_chain_plots.carbayes.diag(chaindat = data[,i], filename = filename[i])
    }
  })
  
}

make_mcmc_chain_plots.carbayes.final<- function(chaindat, filename){
  
    plot(chaindat[[1]][[1]])
    
    jpeg(filename, height = 8, width = 11, units = "in", res=200)
    par(mfrow=c(2,2))
    plot(chaindat[[1]][[1]])
    graphics.off()
  
}

wrap_chain_plotter.carbayes.final <- function(final.dir, chains){
  # this function does not return anything
  # it is internally making plots
  
  purrr::pmap(chains, function(data, name){
    
    filename <- paste0(final.dir, name, "_", colnames(data), ".jpg") 
    
    for(i in 1:ncol(data)){
      make_mcmc_chain_plots.carbayes.final(chaindat = data[,i], filename = filename[i])
    }
  })
  
}




#---------------------------------------------------
# Internal use function, not good for corner cases
#----------------------------------------------------
#' @details Extract betas and covariance parameters from PrevMap.Bayes Object and 
#' calculate the effective sampling size base on my particular mcmc design

get_diag_summ_results.Bayes.PrevMap <- function(mcmc, hpd.coverage = 0.95){
  # coverage
  summresults <- summary(mcmc, hpd.coverage)
  betas <- as.data.frame(summresults$beta )
  betas$param <- rownames(betas)
  covars <- as.data.frame( rbind(summresults$sigma2, summresults$phi, summresults$tau2) )
  covars$param <- c("sigma^2", "phi", "tau^2")
  summresults <- rbind.data.frame(betas, covars) %>% 
    dplyr::select(c("param", dplyr::everything()))

  
  # effective sampling
  params <- summresults$param
  chains <- lapply(params, function(x) return(mcmc$estimate[,x]))
  neff <- lapply(chains, function(x){ return(coda::effectiveSize(coda::mcmc(x)))})
  neff <- data.frame(param = params, neff = unlist(neff))
  
  summresults <- dplyr::left_join(summresults, neff, by = "param")
  
  # trace plots
  plot.new()
  beta <- recordPlot(PrevMap::trace.plot(mcmc, param = "beta", component.beta = 1))
  sigma <- recordPlot(PrevMap::trace.plot(mcmc, param = "sigma2"))
  phi <- recordPlot(PrevMap::trace.plot(mcmc, param = "phi"))
  tau <- recordPlot(PrevMap::trace.plot(mcmc, param = "tau2"))
  traceplots <- list(beta, sigma, phi, tau)
  
  # autocorr plots
  betacorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "beta", component.beta = 1))
  sigmacorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "sigma2"))
  phicorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "phi"))
  taucorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "tau2"))
  autocorr <- list(betacorr, sigmacorr, phicorr, taucorr)
  
  
  ret <- list(summresults = summresults,
              traceplots = traceplots,
              autocorrplots = autocorr)
  return(ret)
  
}



