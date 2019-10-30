library(ggplot2)

#----------------------------------------------------------------------------------------------------
# Basic epi
#----------------------------------------------------------------------------------------------------

logit <- function(x, tol=1e-4){ 
    return( log(((x+tol)/(1-x+tol))) )
}

expit <- function(x, tol=1e-4){ 
  return( 1/(1+exp(-x + tol)) )
}

#----------------------------------------------------------------------------------------------------
# Make final Survey Object to account for DHS Survey Weights
#----------------------------------------------------------------------------------------------------
makecd2013survey <- function(survey = dt){
  # Notes on Lonely PSUs
  # http://r-survey.r-forge.r-project.org/survey/exmample-lonely.html
  options(survey.lonely.psu="adjust")
  dtsrvy <- survey %>% srvyr::as_survey_design(ids = hv001, 
                                               strata = hv023, 
                                               weights = hiv05_wi)
  return(dtsrvy)
}

#----------------------------------------------------------------------------------------------------
# Change return from scale
#----------------------------------------------------------------------------------------------------
my.scale <- function(x, ...){
  ret <- scale(x, ...)
  if(ncol(ret) == 1){
    return(as.numeric(ret)) # coerce matrix of 1 col to numeric vector
  } else{
    stop("Matrix has more than one column")
    }
}

#----------------------------------------------------------------------------------------------------
# Make MLR Task
#----------------------------------------------------------------------------------------------------

make_class_task <- function(data = data, 
                            type = type,
                            target = target,
                            positive = positive,
                            coordinates = NULL){
  
  if(type == "binary"){
    task <- mlr::makeClassifTask(data = data, 
                                 target = target,
                                 positive = positive,
                                 coordinates = coordinates)
  } else if(type == "continuous"){
    task <- mlr::makeRegrTask(data = data, 
                              target = target,
                              coordinates = coordinates)
  }
  
  return(task)
  
}


#----------------------------------------------------------------------------------------------------
# ViVID Epi Theme
#----------------------------------------------------------------------------------------------------

vivid_theme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14), 
                     axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12), 
                     axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11), 
                     legend.position = "bottom",
                     legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
                     legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"),
                     panel.grid = element_blank(),
                     panel.border = element_blank())


