

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









#' @param dag daggity object
#' @param exposure character; 
#' @param outcome character;
#' @param liftoverdf dataframe; a dataframe that matches the daggity names to
#' the actual names in your dataset
#' 
#' @description Purpose of this function is to take the daggity names 
#' and find the matching covariates in the dataset 


get_minimal_set <- function(dag, exposure, outcome, 
                              type = "minimal", effect = "total", 
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
  

  adj.min.count <- min( sapply(adj, length) )
  adj.min <- adj[ sapply(adj, length) == adj.min.count ]
  
  # error catch if more than one minimal set of same lenght
  if(length(adj.min) > 1){
    # randomly pick one
    adj.min <-  adj.min[ sample(1:length(adj.min), size=1) ]
  }
  
  set <- tibble::tibble(dagcovar = as.character(adj.min))
  
  
  # do liftover
  adjset <- dplyr::left_join(set, liftoverdf, by = "dagcovar")
  
  return(unlist(adjset$dhscovar))
  
  
  
}
