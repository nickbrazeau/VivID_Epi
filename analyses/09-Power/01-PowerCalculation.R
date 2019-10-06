library(tidyverse)
library(pwr)
library(furrr)
source("R/00-functions_basic.R")
#...............................................................
# Power Function
#...............................................................
#' @param n numeric; total number of people in population to simulate
#' @param exp_prob numeric; probability of exposure in the population
#' @param p numeric; probability of infection/prevalence of outcome 
#' @param p0 numeric; prevalence among unexposed/probability of outcome among unexposed
powercalculator.glmRR <- function(n=15811, exp_prob=0.5, p=0.03, p0=0.02){
  
  df <- data.frame(obs=factor(seq(1:n)),
                   exp=sample(x=c(0,1), size=n, replace = T, prob=c(exp_prob, 1-exp_prob))) # df of exposure
  p <- 2*p # inv average prev for both groups 
  p0 <- p0 # prev among unexposed
  p1 <- p-p0 # prev among exposed
  RR <- exp(log(p1) - log(p0))

  
    df$dz[df$exp == 1] <- rbinom(sum(df$exp == 1),1,p1)
    df$dz[df$exp == 0] <- rbinom(sum(df$exp == 0),1,p0)

    mod <- glm(dz ~ exp, data=df,
                  family=binomial(link="log"))

    pi <- broom::tidy(mod)$p.value[2]

    ret <- data.frame(RR=RR, p=pi)

    return(ret)

}




#...............................................................
# Make Data Frame for params
#...............................................................
### run lots of these at different levels of p0
p0sim <- seq(0.01, 0.032, by=0.0001)
expprob <- c(0.1, 0.25, 0.5)
exppo <- expand.grid(expprob, p0sim)
poweriters.paramsdf <- tibble::tibble(
  n = 15859, # total pop
  p = 0.03, # prev in population
  exp_prob = exppo[,1],
  p0 = exppo[,2]
) 

# iters to run
iters <- 1e3 
poweriters.paramsdf <- parallel::mclapply(1:iters, function(x) return(poweriters.paramsdf)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(exp_prob, p0)




#...............................................................
# Run in parallel on slurm
#...............................................................
# for slurm on LL
setwd("analyses/09-Power/")
ntry <- 1028 # max number of nodes
sjob <- rslurm::slurm_apply(f = powercalculator.glmRR, 
                            params = poweriters.paramsdf, 
                            jobname = 'powercalc',
                            nodes = ntry, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))

cat("*************************** \n Submitted Power Calc Models \n *************************** ")



