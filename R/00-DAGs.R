

#' @param dag daggity object
#' @param exposure character; 
#' @param outcome character;
#' @param liftoverdf dataframe; a dataframe that matches the daggity names to
#' the actual names in your dataset
#' 
#' @description Purpose of this function is to take the daggity names 
#' and find the matching covariates in the dataset 


get_canonical_set <- function(dag, exposure, outcome, 
                              type = "canonical", effect = "total", 
                              ...,
                              liftoverdf){
  
  # liftover from dhscovar to dagcovar
  exposure <- inner_join(tibble::tibble(dhscovar = exposure), liftoverdf)$dagcovar
  # liftover from dhscovar to dagcovar
  outcome <- inner_join(tibble::tibble(dhscovar = outcome), liftoverdf)$dagcovar
  
  
  adj <- adjustmentSets(dag, 
                        exposure = exposure,
                        outcome = outcome,
                        type = type,
                        effect = effect)
  
  set <- tibble::tibble(dagcovar = as.character(adj[[1]]))
  
  
  # do liftover
  adjset <- dplyr::left_join(set, liftoverdf, by = "dagcovar")
  
  return(unlist(adjset$dhscovar))
  
  
  
}
