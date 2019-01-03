
vivid_map_theme <- theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 13),
                         legend.position = "bottom",
                         legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.85, size = 12),
                         legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, angle = 90, size = 10),
                         axis.text = element_blank(),
                         axis.line = element_blank(),
                         panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent"),
                         panel.grid = element_blank(),
                         panel.border = element_blank())



#----------------------------------------------------------------------------------------------------
#                                           Data Wrangling
#----------------------------------------------------------------------------------------------------

#..................................
# Tidy and DataWrangling functions
#..................................
merge_pr_plsmdm_gemtdt <- function(pr = arpr, plsmdm = panplasmpcrres, ge = ge){
  
  ret <- dplyr::left_join(plsmdm, pr, by="hivrecode_barcode") %>%
    dplyr::left_join(., ge, by = "hv001")
  
  return(ret)
}

#----------------------------------------------------------------------------------------------------
#                                           Spatial Analyses
#----------------------------------------------------------------------------------------------------
#...................................
# Prevalence Summarize, point estimates
#...................................
prev_point_est_summarizer <- function(data, maplvl, plsmdmspec, sfobj){
  
  
  # catch error
  # if( !( any(colnames(x) %in% c("hv001")) & any(colnames(x) %in% c("hv005")) ) ){
  #   stop("Must include cluster and cluster weight. need to check if you are using pr or hiv sample weights")
  # }
  # rlang
  maplvl <- enquo(maplvl)
  plsmdmspec <- enquo(plsmdmspec)
  
  # clusters are weighted (each individual has same weight in cluster)
  ret <- data %>%
    srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>%
    dplyr::group_by(!!maplvl) %>%
    dplyr::summarise(plsmdn = srvyr::survey_total(),
                     plsmd = srvyr::survey_mean(!!plsmdmspec, na.rm = T, vartype = c("se", "ci"), level = 0.95)) %>%
    dplyr::left_join(., sfobj) # attach spatial data, let R figure out the common var
  # return
  return(ret)
}
