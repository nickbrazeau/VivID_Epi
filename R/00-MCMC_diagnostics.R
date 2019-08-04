# TODO add more diagnostics
# https://www2.stat.duke.edu/courses/Fall09/sta290/Lectures/Diagnostics/param-diag.pdf
# consider adding ACF plots
# effectiveSize(theta.MCMC)



# checkConvergence
# calculates Geweke statistic from a series of burn-in and sampling draws. Report whether burn-in length was sufficient based on this statistic.
# (not exported)

checkConvergence <- function(burnin, samples) {
  
  # get number of burnin and sampling iterations
  nburnin <- length(burnin)
  nsamples <- length(samples)
  
  # calculate Geweke diagnostic on combined chain
  chain <- coda::mcmc(c(burnin, samples))
  geweke_z <- coda::geweke.diag(chain, frac1 = nburnin/(nburnin+nsamples), frac2 = nsamples/(nburnin+nsamples))$z
  
  if(is.na(geweke_z)){
    stop("NaN p-value was calculated from Geweke statistic")
  }
  
  # convert to p-value
  geweke_p <- 2*pnorm(abs(geweke_z), lower.tail=FALSE)
  
  # report convergence
  if (geweke_p > 0.05) {
    cat(paste0("convergence reached within defined burn-in period (Geweke p=", round(geweke_p, 3), ")"))
  } else {
    cat(paste0("WARNING: convergence not reached within defined burn-in period (Geweke p=", round(geweke_p,3), ")"))
  }
  
}




#----------------------------------------------
# Internal use function, not good for corner cases
#----------------------------------------------

make_mcmc_chain_plots <- function(chaindat, filename){
  # check if chaindat wasn't fit, i.e. rho in ICAR model
  if(is.na(chaindat[[1]][[1]])){
    plot(1)
  } else {
    
    jpeg(filename, height = 8, width = 11, units = "in", res=200)
    par(mfrow=c(2,2))
    plot(chaindat[[1]][[1]])
    plot(chaindat[[1]][[2]])
    plot(chaindat[[1]][[3]])
    plot(chaindat[[1]][[4]])
    graphics.off()
  }
}

wrap_chain_plotter <- function(tempdir, chains){
  # this function does not return anything
  # it is internally making plots
  
  purrr::pmap(chains, function(data, name){
    
    filename <- paste0(mytempdir, name, "_", colnames(data), ".jpg") 
    
    for(i in 1:ncol(data)){
      make_mcmc_chain_plots(chaindat = data[,i], filename = filename[i])
    }
  })
  
}









