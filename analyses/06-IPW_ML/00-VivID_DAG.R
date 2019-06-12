library(tidyverse)
library(dagitty)
source("R/00-DAGs.R")

vividdag <- dagitty::downloadGraph(x = "dagitty.net/mjKNQhB")

dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))

# iptw sets, going to call these "treatments"
txs <- dcdr %>% 
  dplyr::filter(iptw_model == "y") %>% 
  dplyr::select(-c("risk_factor_raw", "risk_factor_model"))

# find canonical sets
dagliftover <- readxl::read_excel(path = "model_datamaps/dag_dhscovar_liftover.xlsx")
txs$adj_set <- purrr::map(txs$column_name,
                          get_canonical_set,
                          dag = vividdag,
                          outcome = "Pv18s",
                          liftoverdf = dagliftover)



saveRDS(txs, file = "model_datamaps/IPTW_treatments.RDS")





  