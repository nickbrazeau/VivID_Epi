#............................................................
# Gelman-Rubin Diagnostic Carbayes
#...........................................................
#' @title Gelman for the CarBayes Models Random Effect
random_effect_gelman.carbayes <- function(diag_chain1, diag_chain2,
                                          diag_chain3, diag_chain4) {
  # get phis
  gr_phis <- rep(NA, ncol(diag_chain1$phi))
  for (i in 1:length(gr_phis)) {
    c1 <- coda::as.mcmc(diag_chain1$phi[,i])
    c2 <- coda::as.mcmc(diag_chain2$phi[,i])
    c3 <- coda::as.mcmc(diag_chain3$phi[,i])
    c4 <- coda::as.mcmc(diag_chain4$phi[,i])
    mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
    gr_phis[i] <- coda::gelman.diag(mcmclist)$psrf[,1]
  }
  
  # tidy up 
  prov_phis <- data.frame(
    provnum = 1:ncol(diag_chain1$phi),
    gr_phi = gr_phis)
  return(prov_phis)
}


#' @title Gelman for the CarBayes Models Covariates
covars_gelman.carbayes <- function(diag_chain1, diag_chain2,
                                   diag_chain3, diag_chain4) {
  # get betas
  gr_betas <- rep(NA, ncol(diag_chain1$beta))
  for (i in 1:length(gr_betas)) {
    c1 <- coda::as.mcmc(diag_chain1$beta[,i])
    c2 <- coda::as.mcmc(diag_chain2$beta[,i])
    c3 <- coda::as.mcmc(diag_chain3$beta[,i])
    c4 <- coda::as.mcmc(diag_chain4$beta[,i])
    mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
    gr_betas[i] <- coda::gelman.diag(mcmclist)$psrf[,1]
  }
  
  # get tau
  c1 <- coda::as.mcmc(diag_chain1$tau2)
  c2 <- coda::as.mcmc(diag_chain2$tau2)
  c3 <- coda::as.mcmc(diag_chain3$tau2)
  c4 <- coda::as.mcmc(diag_chain4$tau2)
  mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
  gr_tau <- coda::gelman.diag(mcmclist)$psrf[,1]
  
  # get rho
  c1 <- coda::as.mcmc(diag_chain1$rho)
  c2 <- coda::as.mcmc(diag_chain2$rho)
  c3 <- coda::as.mcmc(diag_chain3$rho)
  c4 <- coda::as.mcmc(diag_chain4$rho)
  mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
  gr_rho <- coda::gelman.diag(mcmclist)$psrf[,1]
  
  
  # tidy up 
  prov_covars <- data.frame(
    covars = c(paste0("covars", 1:length(gr_betas)),
               "tau2", "rho"),
    gr_covar = c(gr_betas, gr_tau, gr_rho))
  return(prov_covars)
}

#............................................................
# Gelman-Rubin Diagnostic Prevmap
#...........................................................
#' @title Gelman for the Prevmap Models Random Effect
random_effect_gelman.prevmap <- function(diag_chain1, diag_chain2,
                                         diag_chain3, diag_chain4) {
  # get S matrix
  gr_random_effects <- rep(NA, ncol(diag_chain1$S))
  
  for (i in 1:length(gr_random_effects)) {
    c1 <- coda::as.mcmc(diag_chain1$S[,i])
    c2 <- coda::as.mcmc(diag_chain2$S[,i])
    c3 <- coda::as.mcmc(diag_chain3$S[,i])
    c4 <- coda::as.mcmc(diag_chain4$S[,i])
    mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
    gr_random_effects[i] <- coda::gelman.diag(mcmclist)$psrf[,1]
  }
  
  # tidy up
  cluster_random_effects <- data.frame(
    clstnum = 1:ncol(diag_chain1$S),
    gelmanrubindiag = gr_random_effects
  )
  return(cluster_random_effects)
}


#' @title Gelman for the Prevmap Models Covariates
covariates_gelman.prevmap <- function(diag_chain1, diag_chain2,
                                      diag_chain3, diag_chain4){
  covarsdim <- colnames(diag_chain1$estimate)
  gr_covars <- rep(NA, length(covarsdim))
  for (i in 1:length(covarsdim)) {
    c1 <- coda::as.mcmc(diag_chain1$estimate[,i])
    c2 <- coda::as.mcmc(diag_chain2$estimate[,i])
    c3 <- coda::as.mcmc(diag_chain3$estimate[,i])
    c4 <- coda::as.mcmc(diag_chain4$estimate[,i])
    mcmclist <- coda::as.mcmc.list(list(c1, c2, c3, c4))
    gr_covars[i] <- coda::gelman.diag(mcmclist)$psrf[,1]
  }
  
  # tidy out
  covars_gf <- data.frame(
    covars = covarsdim,
    gelmanrubindiag = gr_covars
  )
  return(covars_gf)
}

#............................................................
# Plotting and Summarizing Tools
#...........................................................
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



