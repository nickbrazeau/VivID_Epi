source("R/00-functions_iptw.R")



covars <- c(
#  "precip_ann_cont_scale_clst",
#     "alt_dem_cont_scale_clst",
#     "wtrdist_cont_scale_clst",
#    "hiv03_fctb",
  "urbanscore_fctb",
  "hlthdist_cont_scale_clst",
#  "anyatm_cont_scale_clst",
#  "pfldh_fctb",
#   "hab57_fctb",
#  "hv104_fctb",
#  "hab1_cont_scale",
#    "hv21345_fctb",
 "wlthrcde_fctb",
 "hv106_fctb"
# "hv246_fctb",
#  "hv014_cont_scale",
#  "hv009_cont_scale",
#  "ITN_fctb" 
    )

outcome <- "ITN_fctb"

# what if i take out altitude
sf::st_geometry(dt) <- NULL

tempdata <- dt %>% 
  dplyr::mutate(urbanscore_fctb = ifelse(urbanscore_cont_scale_clst >= 0, "urban", "rural"),
                urbanscore_fctb = factor(urbanscore_fctb)) %>% 
  dplyr::select(covars, outcome, "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.))

tempdata.dt <- tempdata %>% 
  dplyr::select(-c("longnum", "latnum"))

tempcoords <- tempdata %>% 
  dplyr::select(c("longnum", "latnum"))


task <- mlr::makeClassifTask(data = tempdata.dt, 
                             target = outcome, 
                          coordinates = tempcoords)

task.over <- oversample(task, rate = 50)

nsmp <- sample(1:nrow(tempdata.dt), floor(0.5*nrow(tempdata.dt)))
tempdata.dt.sub <- tempdata.dt[nsmp, ]
tempcoords.sub <- tempcoords[nsmp, ]
task.sub <- mlr::makeClassifTask(data = tempdata.dt.sub, 
                             target = outcome, 
                             coordinates = tempcoords.sub)



# task.over <- oversample(task, rate = 2)
base.learners = c("classif.logreg",
                  # "classif.glmnet", 
                  "classif.kknn",
                  "classif.gamboost",
                  "classif.svm",
                  "classif.randomForest")
baselearners <- lapply(base.learners, makeLearner, predict.type = "prob")
lrn <-  makeStackedLearner(baselearners,
                           predict.type = "prob",
                           method = "hill.climb")





task <- mlr::makeRegrTask(data = tempdata.dt, 
                          target = outcome, 
                          coordinates = tempcoords)

lrn <-  makeStackedLearner(base.learners = c("regr.lm",
                                             #                                          "regr.glmnet", 
                                             "regr.kknn",
                                             "regr.gamboost",
                                             "regr.svm",
                                             "regr.randomForest"),
                           predict.type = "response",
                           method = "average")









mod <- train(lrn, task)
pred <- predict(mod, task)
iptw <- get_iptw_prob(task, pred)
summary(iptw)
hist(iptw)




# map weights

cbind.data.frame(tempdata, iptw) %>% 
  ggplot() + 
  geom_jitter(aes(x=longnum, y=latnum, size=iptw, color = wlthrcde_fctb), alpha = 0.2) 


