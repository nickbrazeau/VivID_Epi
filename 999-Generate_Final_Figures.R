#----------------------------------------------------------------------------------------------------
# Purpose of this script is to generate final figures for Publication
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(rnaturalearth)
library(cowplot)
source("R/00-functions_basic.R")
#......................
# Import Results
#......................
load("data/map_bases/vivid_maps_bases.rda")
load("results/basic_maps_results.rda")
prevmaprasterplots <- readRDS(file = "results/prevmap_raster_plots.rds")



#----------------------------------------------------------------------------------------------------
# MAIN FIGURES
#----------------------------------------------------------------------------------------------------

#----------------------------------------------
# Figure 1B
#----------------------------------------------
pvcasen <- case_n_maps[[2]] +
  ggtitle("") +
  prettybasemap_nodrc + 
  theme(
    legend.position = "right",
    legend.text = element_text(face = "bold", angle = 0, vjust = 0.5, hjust = 0.5)
  )


# make world map
DRC <- readRDS("data/map_bases/gadm/gadm36_COD_0_sp.rds") %>% 
  sf::st_as_sf(.)

bb <- sf::st_bbox(
  sf::st_sf(geom = sf::st_sfc(
    sf::st_point(c(-19, -37)), # left lower
    sf::st_point(c(51, -37)), # right lower
    sf::st_point(c(-19, 38)), # left upper
    sf::st_point(c(51, 38)), # right upper
    crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))
  ))

africa <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>% 
  dplyr::filter(continent == "Africa") %>% 
  sf::st_crop(., y=bb)


africaplot <- ggplot() + 
  geom_sf(data = africa, fill = "#f0f0f0") +
  geom_sf(data = DRC, fill = "#636363")  +
  coord_sf(datum=NA) +
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", size = 1))


# svglite::svglite(file = "results/figures/Figure1B.svg")
jpeg("results/figures/Figure1B.jpg", width = 11, height = 8, units = "in", res = 500)
cowplot::ggdraw() +
  cowplot::draw_plot(pvcasen, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(africaplot, x = 0.02, y= 0.68, width = 0.35, height = 0.25)
graphics.off()

#----------------------------------------------------------------------------------------------------
# SUPPLEMENTAL FIGURES
#----------------------------------------------------------------------------------------------------
#----------------------------------------------
# Bivariate Analysis
#----------------------------------------------
load("results/bivariate_model_results.rda")
covarmap <- readxl::read_excel("model_datamaps/sub_DECODER_covariate_map.xlsx")


pvriskest <- pvivrskfctr_models$glmlog_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::rename(column_name = term) %>% 
  dplyr::mutate(
    column_name = ifelse(grepl("_fctb", column_name),
                         stringr::str_extract(column_name, "[ -~]+_fctb"),
                         column_name)
  )

pvriskest <- dplyr::left_join(pvriskest, covarmap, by = "column_name")
orderrf <- pvriskest %>% 
  dplyr::arrange(level) %>% 
  dplyr::select(abridged_var_label) %>% 
  unlist(.)

pv_bivar_rf_plot <- pvriskest %>% 
  dplyr::mutate(abridged_var_label = factor(abridged_var_label, levels = orderrf, ordered = T)) %>% 
  ggplot() +
  geom_hline(yintercept = 1, color = "#cb181d", linetype = "dashed") +
  geom_pointrange(aes(x = abridged_var_label, y = estimate, 
                      ymin = conf.low, ymax = conf.high,
                      color = factor(level))) +
  scale_color_manual("Level", values = c("#0868ac", "#4eb3d3")) +
  coord_flip() + 
  ggtitle(expression(paste(bold("Bivariate Risk Factor Estimates for "), bolditalic("P. vivax")))) +
  
  ylab("Risk Ratio") + 
  theme(
    plot.title =  element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 14),
    axis.title.x = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 12),
    axis.text = element_text(family = "Helvetica", vjust = 0.5, hjust = 0.5, size = 11),
    axis.title.y = element_blank(),
    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 12),
    legend.text = element_text(family = "Helvetica", vjust = 0.5, hjust = 0.5, size = 10, angle = 0),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank())
    
# PF
pfriskest <- pfalrskfctr_models$glmlog_tidy %>% 
  bind_rows() %>% filter(term != "(Intercept)") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::rename(column_name = term) %>% 
  dplyr::mutate(
    column_name = ifelse(grepl("_fctb", column_name),
                         stringr::str_extract(column_name, "[ -~]+_fctb"),
                         column_name)
  )

pfriskest <- dplyr::left_join(pfriskest, covarmap, by = "column_name")
orderrf <- pfriskest %>% 
  dplyr::arrange(level) %>% 
  dplyr::select(abridged_var_label) %>% 
  unlist(.)

pf_bivar_rf_plot <- pfriskest %>% 
  dplyr::mutate(abridged_var_label = factor(abridged_var_label, levels = orderrf, ordered = T)) %>% 
  ggplot() +
  geom_hline(yintercept = 1, color = "#cb181d", linetype = "dashed") +
  geom_pointrange(aes(x = abridged_var_label, y = estimate, 
                      ymin = conf.low, ymax = conf.high,
                      color = factor(level))) +
  scale_color_manual("Level", values = c("#006d2c", "#41ae76")) +
  coord_flip() + 
  ggtitle(expression(paste(bold("Bivariate Risk Factor Estimates for "), bolditalic("P. falciparum")))) +
  
  ylab("Risk Ratio") + 
  theme(
    plot.title =  element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 14),
    axis.title.x = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 12),
    axis.text = element_text(family = "Helvetica", vjust = 0.5, hjust = 0.5, size = 11),
    axis.title.y = element_blank(),
    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 12),
    legend.text = element_text(family = "Helvetica", vjust = 0.5, hjust = 0.5, size = 10, angle = 0),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank())


# svglite::svglite(file = "results/figures/Figure1B.svg")
jpeg("results/figures/Covar_SuppFig.jpg", width = 11, height = 8, units = "in", res = 500)
cowplot::plot_grid(pv_bivar_rf_plot, 
                   pf_bivar_rf_plot, 
                   nrow = 1)
graphics.off()


#----------------------------------------------------------------------------------------------------
# Additional Figures
#----------------------------------------------------------------------------------------------------

############################################################
#########               MORE MAPS                  ######### 
############################################################
pntestmaps_lrg <- map(pntestmaps, function(x){return(x + prettybasemap_nodrc)})
caseprevmaps_lrg <- map(caseprevmaps, function(x){return(x + prettybasemap_nodrc)})
case_n_maps_lrg <- map(case_n_maps, function(x){return(x + prettybasemap_nodrc)})
prevmaprasterplots_lrg <- map(prevmaprasterplots, function(x){return(x + prettybasemap_nodrc)})


##############################
# Pf, Pv, Po All together
##############################
jpeg("results/figures/pv-pf-po_crude_maps3x3.jpg", width = 11, height = 8, units = "in", res = 500)
gridExtra::grid.arrange(
  pntestmaps_lrg[[1]], caseprevmaps_lrg[[1]], prevmaprasterplots[[1]],
  pntestmaps_lrg[[2]], caseprevmaps_lrg[[2]], prevmaprasterplots[[2]],
  pntestmaps_lrg[[3]], caseprevmaps_lrg[[3]], prevmaprasterplots[[3]],
  ncol=3, top=grid::textGrob("Malaria Species Prevalence in CD2013 DHS", 
                             gp=grid::gpar(fontsize=15, fontfamily = "Arial", fontface = "bold"))) 

graphics.off()

##############################
# Just Pv 
##############################
jpeg("results/figures/pv-crude_maps2x2.jpg", width = 11, height = 8, units = "in", res = 500)
gridExtra::grid.arrange(
  pntestmaps_lrg[[2]], prevmaprasterplots[[2]], 
  caseprevmaps_lrg[[2]], case_n_maps_lrg[[2]],
  ncol=2, top=grid::textGrob("P. vivax Prevalence in CD2013 DHS", 
                             gp=grid::gpar(fontsize=15, fontfamily = "Arial", fontface = "bold"))) 

graphics.off()


##############################
# Ape Overlap
##############################
jpeg("results/figures/pv-ape-overlap-crude_maps1x1.jpg", width = 11, height = 8, units = "in", res = 500)
aperange_nhapv + prettybasemap_nodrc
graphics.off()











