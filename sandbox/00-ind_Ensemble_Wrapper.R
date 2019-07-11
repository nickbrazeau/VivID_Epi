
#' @description 
#' Although a bit unconventional, the IPTW approach does not need to be optimized
#' on "out-of-sample" prediction accuracy. For example, an individual with Pe 
#' may really need to be "down"-weighted becuase of their confounding profile
#' in order to make balanced groups. As a result, the predictive probability and 
#' the balance between the two groups is what we wish to maximize. 
#' 
#' Again, in a bit of an unconventional fashion, I want to stack a few
#' individual learners together instead of using and tuning a stackedlearner
#' as shown here: https://github.com/mlr-org/mlr/issues/1266
#' The main issue is that I cannot perform spatial CV on a stacked object
#' 
#' 
#' I do also want to perform hyper-parameter tuning for 
#' our various learners in order to maximize their indivudal performance.   
#' 
#' This wrapper breaks our ML approach into three processes:
#' First: Set up the resampling strategy. Note, 




# https://github.com/mlr-org/mlr/pull/1041
# https://github.com/mlr-org/mlr/issues/1266
# http://www.cs.uwyo.edu/~larsko/slides/idir18.pdf

make_Ind_Stack <- function(
  
  
)
  #----------------------------------------
  # Resampling Section
  #----------------------------------------

  resamp_sp <- mlr::makeResampleDesc(method = resample$method, 
                                     fold = resample$fold,
                                     reps = resample$reps)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
