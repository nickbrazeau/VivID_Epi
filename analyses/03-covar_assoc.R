
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
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate( risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model) )

# grab risk factors
rskfctr <- dcdr %>% 
  dplyr::filter(risk_factor_model == "y" ) %>% 
  dplyr::select("column_name") %>% 
  unlist(.) %>% 
  unname(.)

dtsrvy <- makecd2013survey(survey = dt)


#------------------------------------------------------------------------------------------
# Analyze for Collinearity 
#------------------------------------------------------------------------------------------
#......................
# Non-parametric Approach
#......................
#######################
# Functions
#######################
my.smd <- function(Var1, Var2, design = dtsrvy){
  
  
  fctrq <- rlang::sym(Var1)

  fctrqdf <- design %>% 
    dplyr::select(!!fctrq) 
  
  if( is.factor(unlist( fctrqdf$variables ) ) ){
    contq <- rlang::sym(Var2)
  } else { # need to reverse now
    contq <- rlang::sym(Var1)
    fctrq <- rlang::sym(Var2)
  }
  
  ret <- design %>% 
    dplyr::mutate(count = 1) %>% 
    dplyr::group_by(!!fctrq) %>% 
    dplyr::summarise(studyn = srvyr::survey_total(count),
                     groupmean = srvyr::survey_mean(!!contq, na.rm = T, vartype = c("var"))
    )
  
  # find min and max if this is a polytomous categorical variable
  ret <- ret %>% 
    dplyr::filter(groupmean %in% c(min(groupmean), max(groupmean)))
  # find smd
  # https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.6607
  ret <- (ret$groupmean[1] - ret$groupmean[2]) / sqrt(  (ret$groupmean_var[1] + ret$groupmean_var[2] ) / 2)
  
  return(ret)
  
}



findcorreval <- function(Var1, Var2){
  var1evalu <- stringr::str_extract(Var1, "cont|fctb|fctm")
  var2evalu <- stringr::str_extract(Var2, "cont|fctb|fctm")
  cb <- cbind(var1evalu, var2evalu)
  cb <- apply(cb, 1, function(x){
    return(
      paste(x[order(x)], collapse = "")
    )
  })

  ret <- unname(
    sapply(cb, function(x){
    switch(x,
           "contcont" = "pearson",
           "fctbfctb" = "chisq",
           "fctmfctm" = "chisq",
           "contfctb" = "smd",
           "contfctm" = "smd",
           "fctbfctm" = "chisq"
           )}) )
  return(ret)
}

docorreval <- function(Var1, Var2, typeeval){
 if(typeeval == "chisq"){
      ret <- survey::svychisq(as.formula(paste0("~", paste(Var1, Var2, sep = "+"))), design = dtsrvy)
      ret <- unname( unlist( ret$statistic ) )
      
 } else if(typeeval == "pearson"){
      ret <- cov2cor( as.matrix( survey::svyvar(as.formula(paste0("~", paste(
        paste0("as.numeric(", Var1, ")"),
        paste0("as.numeric(", Var2, ")"),
        sep = "+"))), 
        design = dtsrvy, na.rm = T) ) )
      # clean up
      ret <- unname( unlist( ret[2,1] ) )
      
 } else if(typeeval == "smd"){
      ret <- my.smd(as.character(Var1), 
                    as.character(Var2), 
                    design = dtsrvy)
      
 } else {
      stop("Error")
 }
  return(ret)
}



decode_corrgrid <- function(decoder = dcdr, longcorrdf){
  decoder1 <- decoder %>% 
    dplyr::select(c("column_name", "var_label")) %>% 
    magrittr::set_colnames(c("Var1", "Var1_label"))
  
  decoder2 <- decoder %>% 
    dplyr::select(c("column_name", "var_label")) %>% 
    magrittr::set_colnames(c("Var2", "Var2_label"))
  
  longcorrdf <- longcorrdf %>% 
    dplyr::left_join(x=., y=decoder1) %>% 
    dplyr::left_join(x=., y=decoder2)
  
  return(longcorrdf)
  
}

#######################
# RUN IT 
#######################

# make corr grid 
corr.grid <- data.frame(t(combn(rskfctr, 2))) %>% 
  magrittr::set_colnames(., c("Var1", "Var2")) %>% 
  dplyr::mutate(
    typeeval = findcorreval(Var1, Var2)
    )
  

corr.grid$corrret <- purrr::pmap(corr.grid, docorreval)
corr.grid$corrret <- unlist(corr.grid$corrret)

corr.grid.decoded <- decode_corrgrid(decoder = dcdr, longcorrdf = corr.grid)

# quick viz for both categoricals
corr.cat.plot <- corr.grid.decoded %>% 
  dplyr::filter(typeeval == "chisq") %>% 
  dplyr::mutate(corrret_scale = my.scale(corrret)) %>% 
  dplyr::select(-c("typeeval", "corrret")) %>% 
  dplyr::mutate(Var1_label = forcats::fct_rev(forcats::fct_reorder(.f = Var1_label, .x = Var1_label, .fun = length)),
                Var2_label = forcats::fct_rev(forcats::fct_reorder(.f = Var2_label, .x = Var2_label, .fun = length))) %>% 
  ggplot() +
  geom_tile(aes(x=Var1_label, y=Var2_label, fill= corrret_scale)) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
  ggtitle("Scaled Chi-Square Statistics for Caterogical Variables") +
  vivid_theme +
  theme(legend.position = "right", 
        legend.text = element_text(angle = 0),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90),
        panel.border = element_blank())


# quick viz for continuous
corr.cont.plot <- corr.grid.decoded %>% 
  dplyr::filter(typeeval == "pearson") %>% 
  dplyr::mutate(corrret_scale = corrret) %>% 
  dplyr::select(-c("typeeval", "corrret")) %>% 
  dplyr::mutate(Var1_label = forcats::fct_rev(forcats::fct_reorder(.f = Var1_label, .x = Var1_label, .fun = length)),
                Var2_label = forcats::fct_rev(forcats::fct_reorder(.f = Var2_label, .x = Var2_label, .fun = length))) %>% 
  ggplot() +
  geom_tile(aes(x=Var1_label, y=Var2_label, fill= corrret_scale)) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
  ggtitle("Pearson Correlations for Continuous Variables") +
  vivid_theme +
  theme(legend.position = "right", 
        legend.text = element_text(angle = 0),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        panel.border = element_blank())


# quick viz for cont-categorical
corr.smd.plot <- corr.grid.decoded %>% 
  dplyr::filter(typeeval == "smd") %>% 
  dplyr::mutate(corrret_scale = my.scale(corrret)) %>% 
  dplyr::select(-c("typeeval", "corrret")) %>% 
  dplyr::mutate(Var1_label = forcats::fct_rev(forcats::fct_reorder(.f = Var1_label, .x = Var1_label, .fun = length)),
                Var2_label = forcats::fct_rev(forcats::fct_reorder(.f = Var2_label, .x = Var2_label, .fun = length))) %>% 
  ggplot() +
  geom_tile(aes(x=Var1, y=Var2, fill= corrret_scale)) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
  ggtitle("Scaled SMD for Continuous-Categorical Variables") +
  vivid_theme +
  theme(legend.position = "right", 
        legend.text = element_text(angle = 0),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        panel.border = element_blank())



plot(corr.cat.plot)
plot(corr.cont.plot)
plot(corr.smd.plot)



#......................
# Parametric Approach
#......................
# VIF
eq <- as.formula(paste0("pv18s~", paste(rskfctr, collapse = "+")))
model.sat <- survey::svyglm(eq,
               design = dtsrvy,
               family = quasibinomial(link="logit"))
summary(model.sat)
vifs <- car::vif(model.sat)
summary(vifs)



#------------------------------------------------------------------------------------------
# Analyze for Missing Data
#-------------------------------------------------------------------------------------------
dtsrvy %>% 
  dplyr::select(rskfctr) %>% 
  as.data.frame(.) %>% 
  mice::md.pattern(., plot=T, rotate.names = T)

# Are my missing cases fundamentally different?
dt.cmpl <- dt %>% 
  dplyr::select(c("hv001", "hv023", "hv005_wi", rskfctr)) %>% 
  dplyr::mutate(
    ismissing = ifelse(rowSums(is.na(.)) == 0, "no", "yes"),
    ismissing = factor(ismissing))

sf::st_geometry(dt.cmpl) <- NULL
dt.cmpl.srvy <- makecd2013survey(dt.cmpl)


misstbl1 <- tableone::svyCreateTableOne(
  data = dt.cmpl.srvy,
  strata = "ismissing",
  vars = rskfctr,
  includeNA = T,
  smd = T,
  test = F)

print(misstbl1, smd = TRUE)

# are missing samples spread out?
dt.cmpl %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::filter(ismissing == "yes") %>% 
  dplyr::summarise(
    n = n()
  ) %>% 
  DT::datatable()


#................................
# Note to any reading this code: 
# This data does not appear to be MCAR... there do appear to be some biases (more or less, poor, rural people with lots of household members)
# However, given that this 0.44% of the unweighted data and 0.43% of the weighted data, we are well within the range of the rule of thumb
# that missing data less than 5-10% (even if it is MNAR) is largely inconsequential (see Bennett 2001). As a result, I feel comfortable subsetting to complete
# cases from here on out without worrying about the effect of any MNAR bias that may or may not be present and would affect my results
#................................

#......................
# Results & Out
#......................
saveRDS(corr.grid, file = "data/derived_data/covariate_correlation_data.rds")








