library(parallelMap)
library(batchtools)

makeClusterFunctionsSlurm(template = "slurm", array.jobs = TRUE,
                          nodename = "localhost", scheduler.latency = 1, fs.latency = 65)


parallelStartBatchtools(bt.resources = list())
parallelLibrary("MASS")
# subsample iris, fit an LDA model and return prediction error
f = function(i) {
  n = nrow(iris)
  train = sample(n, n/2)
  test = setdiff(1:n, train)
  model = lda(Species~., data=iris[train,])
  pred = predict(model, newdata=iris[test,])
  mean(pred$class != iris[test,]$Species)
  Sys.sleep(60)
}

y = parallelMap(f, 1:2)
parallelStop()