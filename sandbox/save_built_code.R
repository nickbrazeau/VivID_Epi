



dt$built_population_2014_logit_cont <- logit(dt$built_population_2014, tol = 1e-3) # tolerance of 1e-3, transform back to real-line with logit 
hist(dt$built_population_2014_logit_cont) # still not normal dist (as expected because of all the 0s) but better

dt <- dt %>% 
  dplyr::mutate(built_population_2014_logit_cont_scale = scale(built_population_2014_logit_cont, center = T, scale = T)) # this becomes (x-mu)/sd
hist(dt$built_population_2014_logit_cont_scale) # still not normal dist (as expected because of all the 0s) but better



dt$nightlights_composite_log_cont_clust <- log(dt$nightlights_composite + 1e-3) # tolerance of 1e-3, not many, many 0s
hist(dt$nightlights_composite_log_cont_clust) # many, many 0s now become tolerance

dt <- dt %>% 
  dplyr::mutate(nightlights_composite_log_cont_clust_scale = scale(nightlights_composite_log_cont_clust, center = T, scale = T)) # this becomes (x-mu)/sd
