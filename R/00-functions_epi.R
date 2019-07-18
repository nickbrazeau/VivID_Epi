#----------------------------------------------------------------------------------------------------
# iptw prob
#----------------------------------------------------------------------------------------------------

get_iptw_prob <- function(task, preds, type){

  if(type == "binary"){
  
   # pull details from mlr for numerator
   pos.class <- mlr::getTaskDesc(task)$positive
   target <- mlr::getTaskTargetNames(task) 
   data <- mlr::getTaskData(task)
   
   
   ps <- mlr::getPredictionProbabilities(preds)
   exposure <- mlr::getPredictionTruth(preds)
   pexp <- mean(exposure == pos.class)
   
   iptw_u <- ifelse(exposure == pos.class,
                    1/ps,
                    1/(1-ps))
   
   iptw_s <- pexp*iptw_u
   
  } else if(type == "continuous"){
    
    preds <- preds$data
    
    # following assumptions in Robbins 2000/Zhu 2015 PMC4749263
    model.num <- lm(preds$truth~1) # this is the intercept in the classic way
    #TODO check if sigma truly necessary 
    # divide residuals by model variance to put on standard normal
    ps.num <- dnorm(
      (preds$truth - model.num$fitted)/summary(model.num)$sigma,
      mean = 0, sd = 1, log = F)
    
    ps.denom <- dnorm(
      (preds$truth - preds$response)/( sd( (preds$truth - preds$response)) ),
      mean = 0, sd = 1, log = F)
    
     iptw_s <- ps.num/ps.denom
     
  } else{
    stop("Type must be either binary or continous")
  }
  
  return(iptw_s)

}



#----------------------------------------------------------------------------------------------------
# glms
#----------------------------------------------------------------------------------------------------

fitsvyglmlogit <- function(outcome, covar){
  dtsrvy_sub <- dtsrvy %>% 
    dplyr::select(c(outcome, covar))
  
  eq <- as.formula(paste0(outcome, "~", covar))
  ret <- survey::svyglm(eq,
                        design = dtsrvy_sub,
                        family = quasibinomial(link="logit"))
  
  return(ret)
  
}

fitsvyglmlog <- function(outcome, covar){
  dtsrvy_sub <- dtsrvy %>% 
    dplyr::select(c(outcome, covar))
  
  eq <- as.formula(paste0(outcome, "~", covar))
  ret <- survey::svyglm(eq,
                        design = dtsrvy_sub,
                        family = quasibinomial(link="log"))
  
  return(ret)
  
}


#----------------------------------------------------------------------------------------------------
# tableone manipulations
#----------------------------------------------------------------------------------------------------

tableone2dataframe <- function(x, columnnames){
  # not table one always sets up neg, pos and then if you have a test as
  # a p-test -- this not corner case proof for t-tests, etc
  capture.output(x <- print(x))
  covarsraw <- factor( c(rownames(x)) ) 
  out <- tibble(
    covars = covarsraw,
    neg = x[,1],
    pos = x[,2])  %>% 
    dplyr::mutate(covars = gsub(" = ", ", ", covars),
                  neg = gsub(" ", "", neg),
                  neg = gsub("\\(", " (", neg),
                  pos = gsub(" ", "", pos),
                  pos = gsub("\\(", " (", pos),
                  matchcol = stringr::str_split_fixed(covarsraw, " ", n = 2)[,1]) %>% 
    magrittr::set_colnames(columnnames)
  
  return(out)
}


mergetableone2table <- function(tableonedf, tabletwoestdf){
  # not built for corner cases obviously

  
  tableonedf <- tableonedf %>% 
    dplyr::mutate(matchcol = ifelse(matchcol == "", NA, matchcol),
                  matchcol = zoo::na.locf(matchcol),
                  lvl = ifelse(
                    grepl("cont|q50", matchcol),
                    "cont",
                    "factor"
                  ),
                  details = stringr::str_split_fixed(tableonedf$Covariates, " ", n = 2)[,2])
  
  # Remember, all continuous variables are scaled in the models but not in the original distributions
  # with the exception of wealth and urbanicity, which are kept as their original scores
  tableonedf <- tableonedf %>% 
    dplyr::mutate(matchcol = gsub("_cont_clst", "_cont_scale_clst", matchcol),
                  matchcol = gsub("_cont$", "_cont_scale", matchcol),
                  matchcol = ifelse(matchcol == "urbanscore_cont_scale_clst", "urbanscore_cont_clst", matchcol),
                  matchcol = ifelse(matchcol == "wlthrcde_combscor_cont_scale", "wlthrcde_combscor_cont", matchcol)
                  )
  
  
  # fix n and wealth from tableone
  tableonedf$lvl[tableonedf$Covariates == "n"] <- "cont"
  
  fctlvls <- stringr::str_extract_all(tableonedf$details[tableonedf$lvl == "factor"], 
                                      "[A-Z]|[a-z]|_", simplify=F)
  fctlvls <- sapply(fctlvls, paste, collapse="") 

  # combine two columns based on ifelse
  tableonedf$matchcol_exp <- NA
  tableonedf$matchcol_exp[tableonedf$lvl == "factor"] <- fctlvls
  tableonedf$matchcol_exp[tableonedf$lvl == "cont"] <- ""
  tableonedf$matchcol <- paste0(tableonedf$matchcol, tableonedf$matchcol_exp)


  # drop extra columns 
  tableonedf <- tableonedf %>% 
    dplyr::select(-c("lvl", "details", "matchcol_exp"))
  
  # bring in table 2
  colnames(tabletwoestdf)[1] <- "matchcol"
  
  ret <- left_join(tableonedf, tabletwoestdf, by = "matchcol")
  
  return(ret)
  }


printriskfactortable2html <- function(rskfcttbl){
  capture.output(x <- print(rskfcttbl))
  knitr::kable(x)
}





