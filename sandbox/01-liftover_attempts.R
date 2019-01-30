#......................
# LIFTOVER
#......................
covar <- readxl::read_excel(path = "internal_datamap_files/liftPR_internal_covar_map.xlsx", sheet = 2)
# dhs strings that need to be NA





test <- dt[, c("hv213", "hv201")]

# start
liftover_covars <- function(covar){
  
  ret <- data.frame(label = attr(covar, "label", exact = TRUE),
                    response = names(attr(covar, "labels", exact = TRUE)),
                    value = unname(attr(covar, "labels", exact = TRUE))) 
  
  return(ret)
  
}

apply(test, 2, liftover_covars)

attr(test$hv213, "label")


dplyr::mutate(recode = ifelse(grepl(pattern = paste(toMatch,collapse="|"), tolower(response)),
                              NA, value))






apply(dt, 2, function(x){stack(attr(x, "labels"))})






t <-  tibble(label = attr(dt$hv201, "label", exact = TRUE),
             response = names(attr(dt$hv201, "labels", exact = TRUE)),
             value = unname(attr(dt$hv201, "labels", exact = TRUE)))  %>% 
  dplyr::mutate(recode = ifelse(grepl(pattern = paste(toMatch,collapse="|"), tolower(response)),
                                NA, value)
  )
t[19,4] <- NA

test$hv201[test$hv201 %in% t$value[is.na(t$recode)]] <- NA





factor(test$hv201, levels = t$recode, labels = t$response[!is.na(t$recode)])




test <- c("missing", "don't know", "don't know brand", "good", "not applicable")
grepl(pattern = paste(toMatch,collapse="|"), test)

apply(dt[,grepl("sh2", colnames(dt))], 2, function(x) sum(!is.na(x))
) 
apply(dt[,grepl("sh23", colnames(dt))], 2, function(x) sum(!is.na(x))
) 
apply(dt[,grepl("hc", colnames(dt))], 2, function(x) sum(!is.na(x))
) 

sum(!is.na(dt$sh240))

t <- if_else(!is.na(dt$ha40), dt$ha40, dt$hb40)