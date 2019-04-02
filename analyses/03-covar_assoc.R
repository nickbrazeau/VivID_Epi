#----------------------------------------------------------------------------------------------------
# Purpose of this script is to examine the associations
# between covariates
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(srvyr)
library(ggcorrplot)
source("R/00-functions_basic.R")


# https://cran.r-project.org/web/packages/jtools/vignettes/svycor.html
# ^ rewrite with weights! 

#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
dcdr <- readxl::read_excel(path = "internal_datamap_files/DERIVED_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))
dtsrvy <- makecd2013survey(survey = dt)


#......................
# Analyze for Corr
#......................
rskfctr <- dcdr$column_name[dcdr$risk_factor_raw == "y"]
dtrskfctr <- dt[,rskfctr]
sf::st_geometry(dtrskfctr) <- NULL
dtrskfctr_exp <- fastDummies::dummy_columns(dtrskfctr, remove_first_dummy = F)
dtrskfctr_exp <- dtrskfctr_exp[, (length(rskfctr)+1):ncol(dtrskfctr_exp)]

# drop one of the binary pairs (they are opposite of each other)
binarykeep <- tibble::tibble(
  column_names = colnames(dtrskfctr_exp),
  base_name = stringr::str_split_fixed(string = column_names, 
                                       pattern = "_(?!.*_)",
                                       n=2)[,1],
  level = stringr::str_split_fixed(string = column_names, 
                                       pattern = "_(?!.*_)",
                                       n=2)[,2]
  ) %>% 
  dplyr::filter(level != "NA") %>% 
  dplyr::filter(grepl("_fctb", base_name)) %>% 
  dplyr::group_by(base_name) %>% 
  dplyr::sample_n(., size = 1, replace = F) 



# sort columns and drop binary 
dtrskfctr_exp <- dtrskfctr_exp %>% 
  dplyr::select(c(
    dplyr::ends_with("NA"),
    dplyr::ends_with("_clst"),
    dplyr::contains("_fctm"),
    binarykeep$column_names
  ))


dtrskfctr.corr <- cor(dtrskfctr_exp)

dtrskfctr.corr_high <- dtrskfctr.corr %>% 
  as.dist(.) %>% 
  broom::tidy(.) %>% 
  dplyr::mutate(abscorr = abs(distance)) %>% 
  dplyr::filter(abscorr >= 0.5 )


corrplot <- ggcorrplot::ggcorrplot(dtrskfctr.corr, 
                       type = "lower",
                       outline.color = "white",
                       colors = c("#313695", "#ffffbf", "#a50026"),
                       hc.order = F) +
  ggtitle("Covariate Correlation Matrix") +
  vivid_theme +
  theme(legend.position = "right", 
        legend.text = element_text(angle = 0),
        axis.text = element_text(size = 5),
        panel.border = element_blank())

#......................
# Look at vif
#......................
rskfctr_ind <- dcdr %>% 
  dplyr::filter(risk_factor_model == "y" & level == "individual") %>% 
  dplyr::filter(! column_name %in% c("insctcd_fctm")) %>%   # takes into account hml20 -- will be perfectly collinear for NO. This is an exploratory var only
  dplyr::select("column_name") %>% 
  unlist(.)

eq <- as.formula(paste0("pv18s~", paste(rskfctr_ind, collapse = "+")))
model.sat <- survey::svyglm(eq,
               design = dtsrvy,
               family = quasibinomial(link="logit"))
summary(model.sat)
car::vif(model.sat)

# old wealth
rskfctr_ind_old_wlth <- rskfctr_ind[!grepl("wlth", rskfctr_ind)]
dt$hv270_fctm <- haven::as_factor(dt$hv270)
xtabs(~dt$hv270_fctm + dt$hv270, addNA = T)

rskfctr_ind_old_wlth <- c(rskfctr_ind_old_wlth, "hv270_fctm")
eq <- as.formula(paste0("pv18s~", paste(rskfctr_ind_old_wlth, collapse = "+")))
dtsrvy <- makecd2013survey(dt)
model.sat_oldwlth <- survey::svyglm(eq,
                            design = dtsrvy,
                            family = quasibinomial(link="logit"))
summary(model.sat_oldwlth)
car::vif(model.sat_oldwlth) # old wealth is even worse



# no house
rskfctr_ind_nohouse <- rskfctr_ind[!grepl("hv21345", rskfctr_ind)]
eq <- as.formula(paste0("pv18s~", paste(rskfctr_ind_nohouse, collapse = "+")))
model.sat.nohouse <- survey::svyglm(eq,
                            design = dtsrvy,
                            family = quasibinomial(link="logit"))
summary(model.sat.nohouse)
car::vif(model.sat.nohouse) # no house nearly halves it

#......................
# Results & Out
#......................
saveRDS(dtrskfctr.corr, file = "data/derived_data/covariate_correlation_data.rds")








