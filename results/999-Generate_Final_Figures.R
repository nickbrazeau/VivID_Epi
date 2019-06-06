#----------------------------------------------------------------------------------------------------
# Purpose of this script is to generate final figures for Publication
#----------------------------------------------------------------------------------------------------
library(tidyverse)

#......................
# Import Results
#......................
load("data/map_bases/vivid_maps_bases.rda")
prevmaprasterplots <- readRDS(file = "results/prevmap_raster_plots.rds")



#----------------------------------------------------------------------------------------------------
# CRUDE MAPS
#----------------------------------------------------------------------------------------------------

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











