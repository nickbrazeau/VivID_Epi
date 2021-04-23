# As desrcibed in the Tusting 2017 manuscript (PMC5319641)
# housing materials are taken into consideration for wealth
# need to recode the wealth variable to avoid controlling for part of our
# effect when considering the covar housing materials (as they did above)
# https://dhsprogram.com/programming/wealth%20index/Steps_to_constructing_the_new_DHS_Wealth_Index.pdf
# note, I have a mix of de jure and de facto but my de jure get subsetted later
#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and recode wealth
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
source("R/00-functions_basic.R")
set.seed(48)

# subset to my study population but include those individuals with missing covariates
# am assuming that the de jure, no geo-located individuals, and hiv weights are a 
# the "true" study population and that the 17 missing covariates are MCAR
dt <- readRDS("data/raw_data/vividpcr_dhs_raw.rds")  %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>%  # drop observations with missing geospatial data 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::filter(hv102 == 1) %>% # subset to de-jure https://dhsprogram.com/data/Guide-to-DHS-Statistics/Analyzing_DHS_Data.htm
  dplyr::filter(hiv05 != 0) # drop observations with samplings weights set to 0
sf::st_geometry(dt) <- NULL

#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hiv05_wi = hiv05/1e6)

# As desrcibed in this manuscript (PMID: 28222094)
# housing materials are  taken into consideration for wealth
# need to recode the wealth variable to avoid controlling for part of our
# effect when considering the covar housing materials (as they did above)
# https://dhsprogram.com/programming/wealth%20index/Steps_to_constructing_the_new_DHS_Wealth_Index.pdf
# PCA to do this 
# Will use the same variables for consideration as (PMID: 28222094)
# (1) source of drinking water; (2) toilet facility; (3) cooking fuel; 
# (4) electricity; ownership of a: (5) radio, (6) television, (7) bicycle, 
# (8) mobile telephone, (9) watch
# all in PR recode

# 1. Drinking Water (categorical: HV201)
# 2. Type of Toilet Facility (HV205) 
# 3. Type of Cooking Fuel (HV226)
# 4. Electricity (HV206)
# 5. Radio (HV207)
# 6. Television (HV208)
# 7. Bicycle (HV210)
# 8. Mobile telephone (HV243A)
# 9. Watch (HV243B)

wlth <- dt %>% 
  as.data.frame(.) %>% # drop tibble and sf classes
  dplyr::select(c(
    "hv201", "hv205", "hv226", 
    "hv206", "hv207", "hv208", "hv210", "hv243a", "hv243b" # binary vars
  )) %>% 
  haven::as_factor(.) 

wlth_fct <- wlth
for(i in 1:ncol(wlth_fct)){
  wlth_fct[,i] <- forcats::fct_drop(wlth_fct[,i])
}

#............................................................
# lift over multifactorial variables
#...........................................................
hv201 <- readxl::read_excel("internal_datamap_files/wlth_recode_hv201.xlsx")
hv205 <- readxl::read_excel("internal_datamap_files/wlth_recode_hv205.xlsx")
hv226 <- readxl::read_excel("internal_datamap_files/wlth_recode_hv226.xlsx")

wlth_fct_binary <- wlth_fct %>% 
  dplyr::left_join(., hv201, by = "hv201") %>% 
  dplyr::left_join(., hv205, by = "hv205") %>% 
  dplyr::left_join(., hv226, by = "hv226") %>% 
  dplyr::select(-c("hv201", "hv205", "hv226"))
summary(wlth_fct_binary)



# RECODE NA values to Missing for Yes-No Original Variables
# This is based on Rutstein point 3b.3
# binary recode
wlth_fct_binary_recode <- wlth_fct_binary
for(i in 1:ncol(wlth_fct_binary_recode)){
  wlth_fct_binary_recode[,i] <- forcats::fct_recode(wlth_fct_binary_recode[,i], no = "missing")
} # some not missing, so will throw warning -- it's ok

# check binary
for(i in 1:ncol(wlth_fct_binary_recode)){
  print(xtabs(~wlth_fct_binary[,i] + wlth_fct_binary_recode[,i], addNA = T))
}

# looks good, no missing/NAs anymore
wlth_fct <- wlth_fct_binary_recode 

# expand factors
wlth_fct_exp <- fastDummies::dummy_columns(wlth_fct)

# sanity
ncol(wlth_fct_exp)
colnames(wlth_fct_exp)

# Tusting's exclusion
# exclude assets where <5 or >95% of households own (per MS)
# 4. Electricity (HV206)
# 5. Radio (HV207)
# 6. Television (HV208)
# 7. Bicycle (HV210)
# 8. Mobile telephone (HV243A)
# 9. Watch (HV243B)
mean(wlth_fct_exp$hv206 == "yes") # ok
mean(wlth_fct_exp$hv207 == "yes") # ok
mean(wlth_fct_exp$hv208 == "yes") # ok
mean(wlth_fct_exp$hv210 == "yes") # ok
mean(wlth_fct_exp$hv243a == "yes") # ok
mean(wlth_fct_exp$hv243b == "yes") # ok

# drop extra columns with factor/chars
wlth_fct_exp <- wlth_fct_exp %>% 
  dplyr::select(c(dplyr::ends_with("_no"), dplyr::ends_with("_yes")))



# add in urbanicity and barcode for split and merge
wlth_fct_exp <- cbind.data.frame(wlth_fct_exp, 
                                 hv025 = haven::as_factor(dt$hv025),
                                 hivrecode_barcode = dt$hivrecode_barcode)

##############################
# RUN  PCA for combined
##############################
# Going to assume we just can use 
# zi1 for the first component, where i is the observations

com1 <- wlth_fct_exp %>% 
  dplyr::select(-c("hv025", "hivrecode_barcode")) %>% 
  prcomp(.)
# compute variance explained
com1$var <- (com1$sdev ^ 2) / sum(com1$sdev ^ 2) * 100
com1$loadings <- abs(com1$rotation)
com1$loadings <- sweep(com1$loadings, 2, colSums(com1$loadings), "/")
# look at var explained
# ~38% in PC1 -- not bad
com1$var
eigen_dir <- cbind.data.frame(nm = rownames(com1$loadings),
                              values = com1$rotation[,1])


#......................
# look at results
#......................
library(ggfortify)
ggplot2::autoplot(com1,
                  loadings = TRUE,
                  loadings.label = TRUE,
                  loadings.label.size  = 2)
# 1. Drinking Water (categorical: HV201) -- lifted over
# 2. Type of Toilet Facility (HV205) -- lifted over
# 3. Type of Cooking Fuel (HV226) -- lifted over
# 4. Electricity (HV206)
# 5. Radio (HV207)
# 6. Television (HV208)
# 7. Bicycle (HV210)
# 8. Mobile telephone (HV243A)
# 9. Watch (HV243B)
# based on the direction of these eigenvectors, this is in *general* lines with my expectations
# own a bicycle is slightly negative but may track with richer adults having a car (or boda boda)

#......................
# tidy into df
#......................
com1df <- data.frame(
  hivrecode_barcode = dt$hivrecode_barcode,
  hv025 = haven::as_factor(dt$hv025),
  com1_scores = com1$x[,1]
)


##############################
# RUN Factor (PCA) analysis for Urban and rural
##############################
urb1 <- wlth_fct_exp %>% 
  dplyr::filter(hv025 == "urban") %>% 
  dplyr::select(-c("hv025", "hivrecode_barcode")) %>% 
  prcomp(.)

urb1df <- cbind.data.frame(
  hivrecode_barcode = dt$hivrecode_barcode[haven::as_factor(dt$hv025) == "urban"],
  urb1_scores = urb1$x[,1]
)

rur1 <- wlth_fct_exp %>% 
  dplyr::filter(hv025 == "rural") %>% 
  dplyr::select(-c("hv025", "hivrecode_barcode")) %>% 
  prcomp(.)

rur1df <- cbind.data.frame(
  hivrecode_barcode = dt$hivrecode_barcode[haven::as_factor(dt$hv025) == "rural"],
  rur1_scores = rur1$x[,1]
)


##############################
# RUN regression and calculate combined score
##############################
wlth_scores <- dplyr::full_join(com1df, rur1df, by = "hivrecode_barcode") %>% 
  dplyr::full_join(., urb1df, by = "hivrecode_barcode")

# check 
xtabs(~wlth_scores$hv025, addNA = T)
sum(!is.na(wlth_scores$rur1_scores))
sum(!is.na(wlth_scores$urb1_scores))

urbmodel <- lm(wlth_scores$com1_scores ~ wlth_scores$urb1_scores)
summary(urbmodel)
urbterms <- broom::tidy(urbmodel)

rurmodel <- lm(wlth_scores$com1_scores ~ wlth_scores$rur1_scores)
summary(rurmodel)
rurterms <- broom::tidy(rurmodel)

# make combined score
wlth_scores <- wlth_scores %>% 
  dplyr::mutate(combscor = 
                  ifelse(hv025 == "rural",
                         rur1_scores*rurterms$estimate[2] + rurterms$estimate[1],
                         ifelse(hv025 == "urban",
                                urb1_scores*urbterms$estimate[2] + urbterms$estimate[1],
                                NA
                         )
                  )
  )
# check
sum(is.na(wlth_scores$combscor))

#...................... 
# look at Combscor, seems reasonable
#...................... 
score_wlth_char <- cbind.data.frame(wlth_fct, data.frame(combscor = wlth_scores$combscor))
summary(score_wlth_char$combscor)
hist(score_wlth_char$combscor)
View(score_wlth_char)
View(eigen_dir)

##############################
# Account for Sampling Weights
##############################
dtsub <- dt %>% 
  dplyr::select(c("hv001", "hv023", "hiv05_wi", "hivrecode_barcode")) %>% 
  dplyr::left_join(., wlth_scores)

dtsrvy <- makecd2013survey(survey = dtsub)
quants <- survey::svyquantile(x = ~combscor, design = dtsrvy, 
                              quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1))

dtsub <- dtsub %>% 
  dplyr::mutate(wlthrcde_fctm = base::cut(x = .$combscor, breaks = quants, 
                                           labels = c("poorest", "poor", 
                                                      "middle", "rich", "richest"),
                                          include.lowest = T)
  ) %>% 
  dplyr::mutate(wlthrcde_fctm = factor(wlthrcde_fctm, levels = c("poorest", "poor", 
                                                                 "middle", "rich", "richest")),
                wlthrcde_fctb = ifelse(wlthrcde_fctm == "poorest", "poor", ifelse(wlthrcde_fctm == "poor", "poor", "not poor")),
                wlthrcde_fctb = factor(wlthrcde_fctb, levels = c("not poor", "poor"))) %>% 
  dplyr::select(c("hivrecode_barcode", "combscor", "wlthrcde_fctm", "wlthrcde_fctb")) %>% 
  dplyr::rename(wlthrcde_combscor_cont = combscor)

summary(dtsub$wlthrcde_combscor_cont)
quants
summary(dtsub$wlthrcde_fctm)
xtabs(~dtsub$wlthrcde_fctm + dtsub$wlthrcde_fctb)

saveRDS(file = "data/derived_data/vividepi_wealth_recoded.rds", object = dtsub)

