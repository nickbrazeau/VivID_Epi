#----------------------------------------------------------------------------------------------------
# Tidy and DataWrangling functions
#----------------------------------------------------------------------------------------------------
merge_pr_plsmdm_gemtdt <- function(pr = arpr, plsmdm = panplasmpcrres, ge = ge){
  
  ret <- dplyr::left_join(plsmdm, pr, by="hivrecode_barcode") %>% 
    dplyr::left_join(., ge, by = "hv001")
  
  return(ret)
}