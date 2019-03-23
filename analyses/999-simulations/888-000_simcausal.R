# https://livefreeordichotomize.com/2019/01/17/understanding-propensity-score-weighting

library(tidyverse)
set.seed(928)
n <- 1000
X <- mvtnorm::rmvnorm(n,
                      mean = c(0.5, 1),
                      sigma = matrix(c(2, 1, 1, 1), ncol = 2)
)

dat <- tibble(
  x_1 = X[, 1],
  x_2 = X[, 2],
  treatment = as.numeric(- 0.5 + 0.25 * x_1 + 0.75 * x_2 + rnorm(n, 0, 1) > 0),
  outcome = as.numeric(rbernoulli(n, plogis(2 * treatment + rnorm(n, 0, 1))))
)

glm(formula = outcome ~ treatment,
    data = dat,
    family = binomial(link = "logit"))

glm(formula = outcome ~ treatment + x_1 + x_2,
    data = dat,
    family = binomial(link = "logit"))