
#............................................................................................................
###### CLASSIFICATION #####
#............................................................................................................
base.learners.classif <- list() 

#....................................
# L1/L2 Regularization 
#....................................
classif.logreg <- mlr::makeLearner("classif.logreg", predict.type = "prob")

classif.lasso <- mlr::makeLearner("classif.cvglmnet", predict.type = "prob")
classif.lasso <- mlr::setHyperPars(classif.lasso, alpha = 1)

classif.elastic <- mlr::makeLearner("classif.cvglmnet", predict.type = "prob")
classif.elastic <- mlr::setHyperPars(classif.elastic, alpha = 0.5)

classif.ridge <- mlr::makeLearner("classif.cvglmnet", predict.type = "prob")
classif.ridge <- mlr::setHyperPars(classif.ridge, alpha = 0)

base.learners.classif <- append(base.learners.classif, list(classif.logreg, 
                                                            classif.lasso, 
                                                            classif.elastic, 
                                                            classif.ridge))

#....................................
# GAM
#....................................
classif.gam <- mlr::makeLearner("classif.gamboost", predict.type = "prob")
base.learners.classif <- append(base.learners.classif, list(classif.gam))
#....................................
# Support Vectors
#....................................
classif.svm.c01g01 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c01g01 <- mlr::setHyperPars(classif.svm.c01g01, cost = 0.01, gamma = 0.01)

classif.svm.c1g01 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c1g01 <- mlr::setHyperPars(classif.svm.c1g01, cost = 1, gamma = 0.01)

classif.svm.c10g01 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c10g01 <- mlr::setHyperPars(classif.svm.c10g01, cost = 10, gamma = 0.01)


classif.svm.c01g1 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c01g1 <- mlr::setHyperPars(classif.svm.c01g1, cost = 0.01, gamma = 1)

classif.svm.c1g1 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c1g1 <- mlr::setHyperPars(classif.svm.c1g1, cost = 1, gamma = 1)

classif.svm.c10g1 <- mlr::makeLearner("classif.svm", predict.type = "prob")
classif.svm.c10g1 <- mlr::setHyperPars(classif.svm.c10g1, cost = 10, gamma = 1)

base.learners.classif <- append(base.learners.classif, list(classif.svm.c01g01, 
                                                            classif.svm.c1g01, 
                                                            classif.svm.c10g01, 
                                                            classif.svm.c01g1, 
                                                            classif.svm.c1g1, 
                                                            classif.svm.c10g1))


#....................................
# Nearest Neighbors
#....................................
classif.kknn5 <- mlr::makeLearner("classif.kknn", predict.type = "prob")
classif.kknn5 <- mlr::setHyperPars(classif.kknn5, k = 5)

classif.kknn10 <- mlr::makeLearner("classif.kknn", predict.type = "prob")
classif.kknn10 <- mlr::setHyperPars(classif.kknn10, k = 10)

classif.kknn25 <- mlr::makeLearner("classif.kknn", predict.type = "prob")
classif.kknn25 <- mlr::setHyperPars(classif.kknn25, k = 25)

classif.kknn50 <- mlr::makeLearner("classif.kknn", predict.type = "prob")
classif.kknn50 <- mlr::setHyperPars(classif.kknn50, k = 50)

base.learners.classif <- append(base.learners.classif, list(classif.kknn5, 
                                                            classif.kknn10, 
                                                            classif.kknn25, 
                                                            classif.kknn50))



#....................................
# Gaussian Process
#....................................
classif.gr.radial <- mlr::makeLearner("classif.gausspr", predict.type = "prob")
classif.gr.radial <- mlr::setHyperPars(classif.gr.radial, kernel = "rbfdot")

classif.gr.laplace <- mlr::makeLearner("classif.gausspr", predict.type = "prob")
classif.gr.laplace <- mlr::setHyperPars(classif.gr.laplace, kernel = "laplacedot")

base.learners.classif <- append(base.learners.classif, list(classif.gr.radial, 
                                                            classif.gr.laplace))

#....................................
# Neural Nets
#....................................
classif.nnet.units1 <- mlr::makeLearner("classif.nnet", predict.type = "prob")
classif.nnet.units1 <- mlr::setHyperPars(classif.nnet.units1, size = 1)

classif.nnet.units2 <- mlr::makeLearner("classif.nnet", predict.type = "prob")
classif.nnet.units2 <- mlr::setHyperPars(classif.nnet.units2, size = 2)

classif.nnet.units3 <- mlr::makeLearner("classif.nnet", predict.type = "prob")
classif.nnet.units3 <- mlr::setHyperPars(classif.nnet.units3, size = 3)

classif.nnet.units4 <- mlr::makeLearner("classif.nnet", predict.type = "prob")
classif.nnet.units4 <- mlr::setHyperPars(classif.nnet.units4, size = 4)

classif.nnet.units5 <- mlr::makeLearner("classif.nnet", predict.type = "prob")
classif.nnet.units5 <- mlr::setHyperPars(classif.nnet.units5, size = 5)

base.learners.classif <- append(base.learners.classif, list(classif.nnet.units1,
                                                            classif.nnet.units2,
                                                            classif.nnet.units3,
                                                            classif.nnet.units4,
                                                            classif.nnet.units5
                                                            ))


#....................................
# Random Forest 
#....................................
classif.rf.m1t500 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m1t500 <- mlr::setHyperPars(classif.rf.m1t500, num.trees  = 500, mtry = 1)                                

classif.rf.m2t500 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m2t500 <- mlr::setHyperPars(classif.rf.m2t500, num.trees  = 500, mtry = 2)    

classif.rf.m3t500 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m3t500 <- mlr::setHyperPars(classif.rf.m3t500, num.trees  = 500, mtry = 3)    

classif.rf.m4t500 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m4t500 <- mlr::setHyperPars(classif.rf.m4t500, num.trees  = 500, mtry = 4)    

classif.rf.m5t500 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m5t500 <- mlr::setHyperPars(classif.rf.m5t500, num.trees  = 500, mtry = 5)    

classif.rf.m1t1000 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m1t1000 <- mlr::setHyperPars(classif.rf.m1t1000, num.trees  = 1000, mtry = 1)                                

classif.rf.m2t1000 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m2t1000 <- mlr::setHyperPars(classif.rf.m2t1000, num.trees  = 1000, mtry = 2)    

classif.rf.m3t1000 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m3t1000 <- mlr::setHyperPars(classif.rf.m3t1000, num.trees  = 1000, mtry = 3)    

classif.rf.m4t1000 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m4t1000 <- mlr::setHyperPars(classif.rf.m4t1000, num.trees  = 1000, mtry = 4)    

classif.rf.m5t1000 <- mlr::makeLearner("classif.ranger", predict.type = "prob")
classif.rf.m5t1000 <- mlr::setHyperPars(classif.rf.m5t1000, num.trees  = 1000, mtry = 5)    



base.learners.classif <- append(base.learners.classif, list(classif.rf.m1t500,
                                                            classif.rf.m2t500,
                                                            classif.rf.m3t500,
                                                            classif.rf.m4t500,
                                                            classif.rf.m5t500,
                                                            
                                                            classif.rf.m1t1000,
                                                            classif.rf.m2t1000,
                                                            classif.rf.m3t1000,
                                                            classif.rf.m4t1000,
                                                            classif.rf.m5t1000
                                                            ))


#............................................................................................................
###### REGRESSION #####
#............................................................................................................
base.learners.regr <- list() 

#....................................
# L1/L2 Regularization 
#....................................
regr.logreg <- mlr::makeLearner("regr.lm", predict.type = "response")

regr.lasso <- mlr::makeLearner("regr.cvglmnet", predict.type = "response")
regr.lasso <- mlr::setHyperPars(regr.lasso, alpha = 1)

regr.elastic <- mlr::makeLearner("regr.cvglmnet", predict.type = "response")
regr.elastic <- mlr::setHyperPars(regr.elastic, alpha = 0.5)

regr.ridge <- mlr::makeLearner("regr.cvglmnet", predict.type = "response")
regr.ridge <- mlr::setHyperPars(regr.ridge, alpha = 0)

base.learners.regr <- append(base.learners.regr, list(regr.logreg, 
                                                      regr.lasso, 
                                                      regr.elastic, 
                                                      regr.ridge))

#....................................
# GAM
#....................................
regr.gam <- mlr::makeLearner("regr.gamboost", predict.type = "response")
base.learners.regr <- append(base.learners.regr, list(regr.gam))
#....................................
# Support Vectors
#....................................
regr.svm.c01g01 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c01g01 <- mlr::setHyperPars(regr.svm.c01g01, cost = 0.01, gamma = 0.01)

regr.svm.c1g01 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c1g01 <- mlr::setHyperPars(regr.svm.c1g01, cost = 1, gamma = 0.01)

regr.svm.c10g01 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c10g01 <- mlr::setHyperPars(regr.svm.c10g01, cost = 10, gamma = 0.01)


regr.svm.c01g1 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c01g1 <- mlr::setHyperPars(regr.svm.c01g1, cost = 0.01, gamma = 1)

regr.svm.c1g1 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c1g1 <- mlr::setHyperPars(regr.svm.c1g1, cost = 1, gamma = 1)

regr.svm.c10g1 <- mlr::makeLearner("regr.svm", predict.type = "response")
regr.svm.c10g1 <- mlr::setHyperPars(regr.svm.c10g1, cost = 10, gamma = 1)

base.learners.regr <- append(base.learners.regr, list(regr.svm.c01g01, 
                                                      regr.svm.c1g01, 
                                                      regr.svm.c10g01, 
                                                      regr.svm.c01g1, 
                                                      regr.svm.c1g1, 
                                                      regr.svm.c10g1))


#....................................
# Nearest Neighbors
#....................................
regr.kknn5 <- mlr::makeLearner("regr.kknn", predict.type = "response")
regr.kknn5 <- mlr::setHyperPars(regr.kknn5, k = 5)

regr.kknn10 <- mlr::makeLearner("regr.kknn", predict.type = "response")
regr.kknn10 <- mlr::setHyperPars(regr.kknn10, k = 10)

regr.kknn25 <- mlr::makeLearner("regr.kknn", predict.type = "response")
regr.kknn25 <- mlr::setHyperPars(regr.kknn25, k = 25)

regr.kknn50 <- mlr::makeLearner("regr.kknn", predict.type = "response")
regr.kknn50 <- mlr::setHyperPars(regr.kknn50, k = 50)

base.learners.regr <- append(base.learners.regr, list(regr.kknn5, 
                                                      regr.kknn10, 
                                                      regr.kknn25, 
                                                      regr.kknn50))



#....................................
# Gaussian Process
#....................................
regr.gr.radial <- mlr::makeLearner("regr.gausspr", predict.type = "response")
regr.gr.radial <- mlr::setHyperPars(regr.gr.radial, kernel = "rbfdot")

regr.gr.laplace <- mlr::makeLearner("regr.gausspr", predict.type = "response")
regr.gr.laplace <- mlr::setHyperPars(regr.gr.laplace, kernel = "laplacedot")

base.learners.regr <- append(base.learners.regr, list(regr.gr.radial, 
                                                      regr.gr.laplace))

#....................................
# Neural Nets
#....................................
regr.nnet.units1 <- mlr::makeLearner("regr.nnet", predict.type = "response")
regr.nnet.units1 <- mlr::setHyperPars(regr.nnet.units1, size = 1)

regr.nnet.units2 <- mlr::makeLearner("regr.nnet", predict.type = "response")
regr.nnet.units2 <- mlr::setHyperPars(regr.nnet.units2, size = 2)

regr.nnet.units3 <- mlr::makeLearner("regr.nnet", predict.type = "response")
regr.nnet.units3 <- mlr::setHyperPars(regr.nnet.units3, size = 3)

regr.nnet.units4 <- mlr::makeLearner("regr.nnet", predict.type = "response")
regr.nnet.units4 <- mlr::setHyperPars(regr.nnet.units4, size = 4)

regr.nnet.units5 <- mlr::makeLearner("regr.nnet", predict.type = "response")
regr.nnet.units5 <- mlr::setHyperPars(regr.nnet.units5, size = 5)

base.learners.regr <- append(base.learners.regr, list(regr.nnet.units1,
                                                      regr.nnet.units2,
                                                      regr.nnet.units3,
                                                      regr.nnet.units4,
                                                      regr.nnet.units5
))


#....................................
# Random Forest 
#....................................
regr.rf.m1t500 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m1t500 <- mlr::setHyperPars(regr.rf.m1t500, num.trees  = 500, mtry = 1)                                

regr.rf.m2t500 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m2t500 <- mlr::setHyperPars(regr.rf.m2t500, num.trees  = 500, mtry = 2)    

regr.rf.m3t500 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m3t500 <- mlr::setHyperPars(regr.rf.m3t500, num.trees  = 500, mtry = 3)    

regr.rf.m4t500 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m4t500 <- mlr::setHyperPars(regr.rf.m4t500, num.trees  = 500, mtry = 4)    

regr.rf.m5t500 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m5t500 <- mlr::setHyperPars(regr.rf.m5t500, num.trees  = 500, mtry = 5)    

regr.rf.m1t1000 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m1t1000 <- mlr::setHyperPars(regr.rf.m1t1000, num.trees  = 1000, mtry = 1)                                

regr.rf.m2t1000 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m2t1000 <- mlr::setHyperPars(regr.rf.m2t1000, num.trees  = 1000, mtry = 2)    

regr.rf.m3t1000 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m3t1000 <- mlr::setHyperPars(regr.rf.m3t1000, num.trees  = 1000, mtry = 3)    

regr.rf.m4t1000 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m4t1000 <- mlr::setHyperPars(regr.rf.m4t1000, num.trees  = 1000, mtry = 4)    

regr.rf.m5t1000 <- mlr::makeLearner("regr.ranger", predict.type = "response")
regr.rf.m5t1000 <- mlr::setHyperPars(regr.rf.m5t1000, num.trees  = 1000, mtry = 5)    



base.learners.regr <- append(base.learners.regr, list(regr.rf.m1t500,
                                                      regr.rf.m2t500,
                                                      regr.rf.m3t500,
                                                      regr.rf.m4t500,
                                                      regr.rf.m5t500,
                                                      
                                                      regr.rf.m1t1000,
                                                      regr.rf.m2t1000,
                                                      regr.rf.m3t1000,
                                                      regr.rf.m4t1000,
                                                      regr.rf.m5t1000
))

remove <- ls()
remove <- remove[! remove %in% c("base.learners.regr", "base.learners.classif")]
remove <- remove[grepl("regr.|classif.", remove)]
rm(list = c(remove, "remove"))
