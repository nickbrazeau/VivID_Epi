
#----------------------------------------------------------------------------------------------------
# Basic epi
#----------------------------------------------------------------------------------------------------

logit <- function(x, tol=1e-4){ 
    return( log(((x+tol)/(1-x+tol))) )
}

#----------------------------------------------------------------------------------------------------
# Make final Survey Object to account for DHS Survey Weights
#----------------------------------------------------------------------------------------------------
makecd2013survey <- function(survey = dt){
  options(survey.lonely.psu="certainty")
  dtsrvy <- survey %>% srvyr::as_survey_design(ids = hv001, 
                                           strata = hv023, 
                                           weights = hv005_wi)
  return(dtsrvy)
}


#----------------------------------------------------------------------------------------------------
# ViVID Epi Theme
#----------------------------------------------------------------------------------------------------

vivid_theme <- theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 14), 
                     axis.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 12), 
                     axis.text = element_text(famil = "Arial", hjust = 0.5, size = 11), 
                     legend.position = "bottom",
                     legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
                     legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent"),
                     panel.grid = element_blank(),
                     panel.border = element_blank())


