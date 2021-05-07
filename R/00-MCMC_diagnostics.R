
#---------------------------------------------------
# Internal use function, not good for corner cases
#----------------------------------------------------

make_mcmc_chain_plots.carbayes <- function(diag.dir, name, chainnum,
                                           betachains, phichains, tau2chains, rhochains){
  
  filename <- paste0(diag.dir, name, "_chain_", chainnum, ".pdf")
  pdf(filename, height = 8, width = 11)
  plot(betachains)
  plot(phichains)
  plot(tau2chains)
  plot(rhochains)
  graphics.off()
  
}



#---------------------------------------------------
# Internal use function, not good for corner cases
#----------------------------------------------------
#' @details Extract betas and covariance parameters from PrevMap.Bayes Object and 
#' calculate the effective sampling size base on my particular mcmc design

get_diag_summ_results.Bayes.PrevMap <- function(mcmc, name, hpd.coverage = 0.95){
  if (grepl("covars", name)) {
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
    # bring together
    summresults <- dplyr::left_join(summresults, neff, by = "param")
    
    # Geweke diagnositc
    geweke.diag <- lapply(chains, function(x){ return(coda::geweke.diag(coda::mcmc(x)))})
    geweke <- data.frame(param = params, Geweke.diag = purrr::map_dbl(geweke.diag, "z"))
    geweke <- geweke %>% 
      dplyr::mutate(Geweke.pval = 2*pnorm(abs(Geweke.diag), lower.tail=FALSE))
    # bring together
    summresults <- dplyr::left_join(summresults, geweke, by = "param")
    
    
    # trace plots
    plot.new()
    beta1 <- recordPlot(PrevMap::trace.plot(mcmc, param = "beta", component.beta = 1))
    beta2 <- recordPlot(PrevMap::trace.plot(mcmc, param = "beta", component.beta = 2))
    beta3 <- recordPlot(PrevMap::trace.plot(mcmc, param = "beta", component.beta = 3))
    sigma <- recordPlot(PrevMap::trace.plot(mcmc, param = "sigma2"))
    phi <- recordPlot(PrevMap::trace.plot(mcmc, param = "phi"))
    tau <- recordPlot(PrevMap::trace.plot(mcmc, param = "tau2"))
    traceplots <- list(beta1, beta2, beta3,
                       sigma, phi, tau)
    dev.off()
    
    plot.new()
    plot(x = 1)
    # autocorr plots
    betacorr1 <- recordPlot(PrevMap::autocor.plot(mcmc, param = "beta", component.beta = 1))
    betacorr2 <- recordPlot(PrevMap::autocor.plot(mcmc, param = "beta", component.beta = 2))
    betacorr3 <- recordPlot(PrevMap::autocor.plot(mcmc, param = "beta", component.beta = 3))
    sigmacorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "sigma2"))
    phicorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "phi"))
    taucorr <- recordPlot(PrevMap::autocor.plot(mcmc, param = "tau2"))
    autocorr <- list(betacorr1, betacorr2, betacorr3,
                     sigmacorr, phicorr, taucorr)
    dev.off()
    
    ret <- list(summresults = summresults,
                traceplots = traceplots,
                autocorrplots = autocorr)
    
  } else if (grepl("intercept", name)) {
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
    # bring together
    summresults <- dplyr::left_join(summresults, neff, by = "param")
    
    # Geweke diagnositc
    geweke.diag <- lapply(chains, function(x){ return(coda::geweke.diag(coda::mcmc(x)))})
    geweke <- data.frame(param = params, Geweke.diag = purrr::map_dbl(geweke.diag, "z"))
    geweke <- geweke %>% 
      dplyr::mutate(Geweke.pval = 2*pnorm(abs(Geweke.diag), lower.tail=FALSE))
    # bring together
    summresults <- dplyr::left_join(summresults, geweke, by = "param")
    
    
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
  } 
  return(ret)
  
}



