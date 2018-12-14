
recode <- tolower(gsub("CD", "", stringr::str_extract(list.files("datasets/"), "CD[A-Z][A-Z]")))

recode = list.files("datasets", full.names = T)[1]

findrecodemtdt <- function(recode){
  df <- readRDS(file = recode)
  name <- gsub("datasets/CD", "", stringr::str_extract(recode, "datasets/CD[A-Z][A-Z]"))
  mtdt <- as.tibble(rdhs::get_variable_labels(df))
  readr::write_csv(mtdt, path = paste0("internal_datamap_files/", "mtdt", name, "_covar_names_labels.csv"))
  
  return(0)
  
}

recodes <- list.files("datasets", full.names = T)[!grepl("GC|GE", recodes)]
lapply(recodes, findrecodemtdt)


