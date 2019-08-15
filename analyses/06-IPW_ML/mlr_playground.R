source("R/00-functions_iptw.R")

<<<<<<< Updated upstream
=======
covars <- getTaskFeatureNames(txs$task[[2]])


covars <- c(
  # "urbanscore_cont_scale_clst", 
   #          "wtrdist_cont_scale_clst",
            "precip_ann_cont_scale_clst"
   #          "alt_dem_cont_scale_clst"
             )
# 
# covars <- covars[!covars %in% c("built_population_2014_cont_scale_clst", "nightlights_composite_cont_scale_clst", "all_population_count_2015_cont_scale_clst", "travel_times_2015_cont_scale_clst" )]
# covars <- covars[grepl("_clst", covars)]
# covars <- c(covars, "urbanscore_cont_scale_clst")
>>>>>>> Stashed changes


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
<<<<<<< Updated upstream
  dplyr::mutate(urbanscore_fctb = ifelse(urbanscore_cont_scale_clst >= 0, "urban", "rural"),
                urbanscore_fctb = factor(urbanscore_fctb)) %>% 
  dplyr::select(covars, outcome, "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.))
=======
  dplyr::mutate(temp_ann_bi_clst = ifelse(temp_ann_cont_clst > 20, "high", "low"),
                temp_ann_bi_clst = factor(temp_ann_bi_clst)) %>% 
  dplyr::select(covars, "temp_ann_bi_clst", "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum"))
>>>>>>> Stashed changes

tempdata.dt <- tempdata %>% 
  dplyr::select(-c("longnum", "latnum"))

<<<<<<< Updated upstream
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
=======
tempcoords <- dt %>% 
  dplyr::select(covars, "temp_ann_cont_clst", "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("longnum", "latnum"))


task <- mlr::makeClassifTask(data = tempdata, target = "temp_ann_bi_clst", 
                          coordinates = tempcoords)
lrn <-  makeStackedLearner(base.learners = c("classif.logreg",
   #                                          "classif.glmnet", 
                                             "classif.kknn",
                                              "classif.gamboost",
                                             "classif.svm",
                                             "classif.randomForest"),
>>>>>>> Stashed changes
                           predict.type = "response",
                           method = "average")









mod <- train(lrn, task)
pred <- predict(mod, task)
iptw <- get_iptw_prob(task, pred)
summary(iptw)
hist(iptw)




# map weights

<<<<<<< Updated upstream
cbind.data.frame(tempdata, iptw) %>% 
  ggplot() + 
  geom_jitter(aes(x=longnum, y=latnum, size=iptw, color = wlthrcde_fctb), alpha = 0.2) 
=======
cbind.data.frame(dt.ml, iptw) %>% 
  ggplot() + 
  geom_point(aes(x=longnum, y=latnum, color = temp_ann_cont_clst ))
  geom_point(aes(x=longnum, y=latnum, color=iptw), alpha = 0.2) 
>>>>>>> Stashed changes


