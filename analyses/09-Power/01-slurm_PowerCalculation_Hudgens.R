library(tidyverse)
library(pwr)
source("R/00-functions_basic.R")
#...............................................................
# Power Function
#...............................................................
#' @param n numeric; total number of people in population to simulate
#' @param exp_prob numeric; probability of exposure in the population
#' @param p numeric; probability of infection/prevalence of outcome 
#' @param p0 numeric; prevalence among unexposed/probability of outcome among unexposed
powercalculator.glmOR <- function(n=15879, exp_prob=0.5, p=0.03, p0=0.02){
  
  df <- data.frame(obs=factor(seq(1:n)),
                   exp=sample(x=c(0,1), size=n, replace = T, prob=c(exp_prob, 1-exp_prob))) # df of exposure
  p <- 2*p # inv average prev for both groups 
  p0 <- p0 # prev among unexposed
  p1 <- p-p0 # prev among exposed
  OR <- exp(logit(p1) - logit(p0))

  
    df$dz[df$exp == 1] <- rbinom(sum(df$exp == 1),1,p1)
    df$dz[df$exp == 0] <- rbinom(sum(df$exp == 0),1,p0)

    mod <- glm(dz ~ exp, data=df,
                  family=binomial(link="logit"))

    pi <- broom::tidy(mod)$p.value[2]

    ret <- data.frame(OR=OR, p=pi)

    return(ret)

}


#...............................................................
# Make Data Frame for params
#...............................................................
### run lots of these at different levels of p0
p0sim <- seq(0.01, 0.032, by=0.0001)
expprob <- c(0.1, 0.25, 0.5)
paramsdf <- tibble::tibble(
  n = 15879, # total pop
  p = 0.03, # prev in population
  exp_prob = rep(expprob, length(p0sim)),
  p0 = rep(p0sim, length(expprob))
) 

# iters to run
iters <- 1e3
paramsdf <- lapply(1:iters, function(x) return(paramsdf)) %>% 
  dplyr::bind_rows() %>% 
  dplyr::arrange(exp_prob, p0)




#...............................................................
# Run the Power Calc on slurm
#...............................................................

# for slurm on LL
setwd("analyses/09-Power/")
ntry <- nrow(paramsdf)

sjob <- rslurm::slurm_apply(f = powercalculator.glmOR, 
                            params = paramsdf, 
                            jobname = 'powercalcs',
                            nodes = 128, 
                            cpus_per_node = 1, 
                            submit = T,
                            slurm_options = list(mem = 16000,
                                                 array = sprintf("0-%d%%%d", 
                                                                 ntry, 
                                                                 128),
                                                 'cpus-per-task' = 8,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))



