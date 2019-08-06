source("R/00-functions_iptw.R")

covars <- getTaskFeatureNames(params$task[[3]])

covars <- c("wtrdist_cont_scale_clst",
            "temp_ann_cont_scale_clst",
            "precip_ann_cont_scale_clst",
            "hv106_fctb",
            "wlthrcde_fctb")

covars <- c("urbanscore_cont_scale_clst")

covars <- covars[!covars %in% c("built_population_2014_cont_scale_clst", "nightlights_composite_cont_scale_clst", "all_population_count_2015_cont_scale_clst", "travel_times_2015_cont_scale_clst" )]
covars <- covars[grepl("_clst", covars)]
covars <- c(covars, "urbanscore_cont_scale_clst")
sf::st_geometry(dt) <- NULL
tempdata <- dt %>% 
  dplyr::select(covars, "hlthdist_cont_clst", "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum"))


tempcoords <- dt %>% 
  dplyr::select(covars, "hlthdist_cont_clst", "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("longnum", "latnum"))


task <- mlr::makeRegrTask(data = tempdata, target = "hlthdist_cont_clst", coordinates = tempcoords)
lrn <-  makeStackedLearner(base.learners = c("regr.lm",
                                             "regr.glmnet", 
                                             "regr.kknn",
                                             "regr.svm",
                                             "regr.randomForest"),
                           predict.type = "response",
                           method = "average")

mod <- train(lrn, task)
pred <- predict(mod, task)
iptw <- get_iptw_prob(task, pred)
summary(iptw)
