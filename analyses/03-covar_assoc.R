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
dcdr <- readxl::read_excel(path = "internal_datamap_files/DECODER_covariate_map.xlsx", sheet = 1) %>% 
  dplyr::mutate( colllinear_level = ifelse(is.na(colllinear_level), "n", colllinear_level) )

# make id framework 
rskfctr.id <- dcdr %>% 
  dplyr::filter(colllinear_level == "id" ) %>% 
  dplyr::select("column_name") %>% 
  unlist(.)

dtsrvy <- makecd2013survey(survey = dt)

# make cluster-level framework
rskfctr.sp <- dcdr %>% 
  dplyr::filter(colllinear_level == "sp" ) %>% 
  dplyr::select("column_name") %>% 
  unlist(.)

clst <- dt[,c("hv001", rskfctr.sp)] %>% 
  dplyr::filter(!duplicated(.))


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
corr.cat.plot <- corr.grid.decoded %>% 
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
  mice::md.pattern(., plot=T)



#......................
# Results & Out
#......................
saveRDS(corr.grid, file = "data/derived_data/covariate_correlation_data.rds")








