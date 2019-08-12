source("R/00-functions_iptw.R")

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




# what if i take out altitude
sf::st_geometry(dt) <- NULL

tempdata <- dt %>% 
  dplyr::mutate(temp_ann_bi_clst = ifelse(temp_ann_cont_clst > 20, "high", "low"),
                temp_ann_bi_clst = factor(temp_ann_bi_clst)) %>% 
  dplyr::select(covars, "temp_ann_bi_clst", "longnum", "latnum") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(-c("longnum", "latnum"))


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
                           predict.type = "response",
                           method = "average")

mod <- train(lrn, task)
pred <- predict(mod, task)
iptw <- get_iptw_prob(task, pred)
summary(iptw)
hist(iptw)




# map weights

cbind.data.frame(dt.ml, iptw) %>% 
  ggplot() + 
  geom_point(aes(x=longnum, y=latnum, color = temp_ann_cont_clst ))
  geom_point(aes(x=longnum, y=latnum, color=iptw), alpha = 0.2) 


