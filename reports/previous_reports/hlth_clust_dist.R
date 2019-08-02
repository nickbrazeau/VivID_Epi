library(tidyverse)
library(mlr)


dt <- readRDS("data/derived_data/vividepi_recode.rds")
sf::st_geometry(dt) <- NULL
mp <- dt %>% 
  dplyr::select(c("hv001", "wtrdist_cont_clst", "wtrdist_cont_scale_clst", "hlthdist_cont_clst", "hlthdist_cont_clst")) %>% 
  dplyr::filter(!duplicated(.))


plot(dt$proximity_to_water ~ dt$wtrdist_cont_clst)
hist(dt$proximity_to_water)
hist(dt$wtrdist_cont_clst)


temp <- lm(hlthdist_cont_clst ~ urbanscore_cont_clst,
           data = dt)
summary(temp)
mean(summary(temp)$residuals^2)


tempdf <- dt %>% 
#  dplyr::select(c("hlthdist_cont_clst", "urbanscore_cont_clst", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "nightlights_composite_scale", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "all_population_count_2015_scale", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "travel_times_2015_scale", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "built_population_2014_scale", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "nightlights_composite", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "all_population_count_2015", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "travel_times_2015", "ITN_fctb", "precip_ann_cont_clst"))
#  dplyr::select(c("hlthdist_cont_clst", "built_population_2014", "ITN_fctb", "precip_ann_cont_clst"))
#
#  dplyr::select(c("hlthdist_cont_clst", "built_population_2014", "all_population_count_2015", "nightlights_composite", "travel_times_2015"))
  dplyr::select(c("hlthdist_cont_clst", "built_population_2014_scale", "all_population_count_2015_scale", "nightlights_composite_scale", "travel_times_2015_scale"))

task = makeRegrTask(data = tempdf, target = "hlthdist_cont_clst")
mod = train(makeLearner("regr.glmnet"), task)
task.pred = predict(mod, task = task)
performance(task.pred)
wi <- get_iptw_prob(task, task.pred) 
hist( wi )


# dats
colnames(task.pred$data)
plot(task.pred$data$truth ~ task.pred$data$response)

mp.left <- mp %>% 
  dplyr::select(c("hv001", "hlthdist_cont_clst"))


dat <- task.pred$data 
dat <- cbind(dat, wi)
colnames(dat)[2] <- "hlthdist_cont_clst"
dat <- left_join(dat, mp.left)

temp <- left_join(tempdf, dat) %>% 
  dplyr::mutate(diff = abs(hlthdist_cont_clst - response)) %>% 
  dplyr::filter(!duplicated(.))


