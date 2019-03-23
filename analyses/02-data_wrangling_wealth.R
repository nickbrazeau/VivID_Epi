#----------------------------------------------------------------------------------------------------
# Purpose of this script is to wrangle and recode wealth
#----------------------------------------------------------------------------------------------------
# libraries and imports
library(tidyverse)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
tol <- 1e-3
set.seed(42)

dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/raw_data/vividpcr_dhs_raw.rds")
#.............
# weights
#.............
dt <- dt %>% 
  dplyr::mutate(hv005_wi = hv005/1e6
  )

# As desrcibed in this manuscript (PMID: 28222094)
# housing materials are  taken into consideration for wealth
# need to recode the wealth variable to avoid controlling for part of our
# effect when considering the covar housing materials (as they did above)
# https://dhsprogram.com/programming/wealth%20index/Steps_to_constructing_the_new_DHS_Wealth_Index.pdf
# Factor analysis/PCA to do this 
# Will use the same variables for consideration as (PMID: 28222094)
# (1) source of drinking water; (2) toilet facility; (3) cooking fuel; 
# (4) electricity; ownership of a: (5) radio, (6) television, (7) bicycle, 
# (8) mobile telephone, (9) watch
# all in PR recode

# 1. Drinking Water (categorical: HV201)
# 2. Type of Toilet Facility (HV205) & Shared Toilet (HV225) from Rutstein point 3b.2.b
# 3. Type of Cooking Fuel (HV226)
# 4. Electricity (HV206)
# 5. Radio (HV207)
# 6. Television (HV208)
# 7. Bicycle (HV210)
# 8. Mobile telephone (HV243A)
# 9. Watch (HV243B)

wlth <- dt %>% 
  dplyr::select(c(
    "hv201", "hv205", "hv226", 
    "hv206", "hv207", "hv208",  "hv225", "hv210", "hv243a", "hv243b" # binary vars
    
  )) %>% 
  haven::as_factor(.) %>% 
  as.data.frame(.) # drop tibble classes

wlth_fct <- wlth
for(i in 1:ncol(wlth_fct)){
  wlth_fct[,i] <- forcats::fct_drop(wlth_fct[,i])
}

# RECODE NA values to Missing for Yes-No Original Variables
# This is based on Rutstein point 3b.3
wlth_fct_binary <- wlth_fct[, c("hv206", "hv207", "hv208",  "hv225", "hv210", "hv243a", "hv243b")]
wlth_fct_multi <- wlth_fct[, c("hv201", "hv205", "hv226")]
# note, hv225 has NAs in it that aren't "missing" going to recode thme to missing
wlth_fct_binary$hv225[is.na(wlth_fct_binary$hv225)] <- factor("missing")

# binary recode
wlth_fct_binary_recode <- wlth_fct_binary
for(i in 1:ncol(wlth_fct_binary_recode)){
  wlth_fct_binary_recode[,i] <- forcats::fct_recode(wlth_fct_binary_recode[,i], no = "missing")
} # some not missing, so will throw warning -- it's ok

# check binary
for(i in 1:ncol(wlth_fct_binary_recode)){
  print(xtabs(~wlth_fct_binary[,i] + wlth_fct_binary_recode[,i], addNA = T))
}

# multinomial recode
wlth_fct_multi_recode <- wlth_fct_multi
# Considering Rutstein point 3b.2.a, it appears these are already recoded together
for(i in 1:ncol(wlth_fct_multi_recode)){
  wlth_fct_multi_recode[,i] <- forcats::fct_recode(wlth_fct_multi_recode[,i], NULL = "missing")
} 
# check multi
for(i in 1:ncol(wlth_fct_multi_recode)){
  print(xtabs(~wlth_fct_multi[,i] + wlth_fct_multi_recode[,i], addNA = T))
}

# RECODE Toilet Facilities
# This is based on Rutstein point 3b.3
wlth_fct <- dplyr::bind_cols(wlth_fct_multi_recode, wlth_fct_binary_recode)
wlth_fct <- wlth_fct %>% 
  dplyr::mutate(toiletfacil = paste(hv205, hv225, sep = "_"),
                toiletfacil = factor( ifelse(grepl("NA", toiletfacil), NA, toiletfacil)) )
xtabs(~wlth_fct$hv205 + wlth_fct$toiletfacil, addNA = T) # have to take into consideration two lost to toilet type missing
wlth_fct <- wlth_fct %>% 
  dplyr::select(-c("hv205", "hv225"))

# impute missing values by "average" -- going to flip a weighted coin/die
# This is based on Rutstein point 4b
impute_missingness <- function(fct){
  if(sum(is.na(fct)) == 0){
    return(fct)
  } else{
    probs <- as.numeric(table(fct))/sum(as.numeric(table(fct))) # purposefully don't include NAs
    options <- names(table(fct))
    cnt <- sum(is.na(fct))
    fct[is.na(fct)] <- sample(x = options,
                              size = cnt,
                              replace = T,
                              prob = probs)
    return(fct)
    
  }
}

wlth_fct <- apply(wlth_fct, 2, impute_missingness)
sum(is.na(wlth_fct)) # good to go

# expand factors
wlth_fct_exp <- fastDummies::dummy_columns(wlth_fct)
# missing observations get coded to all 0s -- ok
strtcol <- min( which(grepl("hv201_", colnames(wlth_fct_exp))) )
wlth_fct_exp <- wlth_fct_exp[, strtcol:ncol(wlth_fct_exp)] # drop original codes

# exclude assets where <55 or >95% of households own (per MS)
dropassets <- wlth_fct_exp %>% 
  tidyr::gather(., key = "assets", value = "owns") %>% 
  dplyr::group_by(assets) %>% 
  dplyr::summarise(pctowns = mean(owns)) %>% 
  dplyr::filter(pctowns < 0.05 | pctowns > 0.95)
# drop
wlth_fct_exp <- wlth_fct_exp[, !( colnames(wlth_fct_exp) %in% dropassets)]


# add in urbanicity and barcode for split and merge
wlth_fct_exp <- cbind.data.frame(wlth_fct_exp, 
                                 hv025 = haven::as_factor(dt$hv025),
                                 hivrecode_barcode = dt$hivrecode_barcode)

##############################
# RUN Factor (PCA) for combined
##############################
# Per this line from Rutstein point 4b
# "principal components extraction using correlation method with one factor extracted"
# am goingt to assume PC1 is extracted and we are calcuating
# zi1 for the first component, where i is the observations

com1 <- wlth_fct_exp %>% 
  dplyr::select(-c("hv025", "hivrecode_barcode")) %>% 
  prcomp(.)
# compute variance explained
com1$var <- (com1$sdev ^ 2) / sum(com1$sdev ^ 2) * 100
# ~28% in PC1 -- not bad

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

##############################
# Account for Sampling Weights
##############################
dtsub <- dt %>% 
  dplyr::select(c("hv001", "hv023", "hv005_wi", "hivrecode_barcode")) %>% 
  dplyr::left_join(., wlth_scores)

dtsrvy <- makecd2013survey(survey = dtsub)
quants <- survey::svyquantile(x = ~combscor, design = dtsrvy, quantiles = c(0, 0.2, 0.4, 0.6, 0.8, 1))

dtsub <- dtsub %>% 
  dplyr::mutate(wlthrcde_fctm = base::cut(x = dtsub$combscor, breaks = quants, 
                                           labels = c("poorest", "poor", "middle", "rich", "richest"))
  ) %>% 
  dplyr::select(c("hivrecode_barcode", "wlthrcde_fctm"))


saveRDS(file = "data/derived_data/vividepi_wealth_recoded.rds", object = dtsub)


