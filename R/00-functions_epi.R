
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
                    grepl("SD|q50|n", matchcol),
                    "cont",
                    "factor"
                  ),
                  details = stringr::str_split_fixed(tableonedf$Covariates, " ", n = 2)[,2])
  
  fctlvls <- stringr::str_extract_all(tableonedf$details[tableonedf$lvl == "factor"], 
                                      "[A-Z]|[a-z]|_", simplify=F)
  fctlvls <- sapply(fctlvls, paste, collapse="") 

  # combine two columns based on ifelse
  tableonedf$matchcol_exp <- NA
  tableonedf$matchcol_exp[tableonedf$lvl == "factor"] <- fctlvls
  tableonedf$matchcol_exp[tableonedf$lvl == "cont"] <- ""
  tableonedf$matchcol = paste0(tableonedf$matchcol, tableonedf$matchcol_exp)
  
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

