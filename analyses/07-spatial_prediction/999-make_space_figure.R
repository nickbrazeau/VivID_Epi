## .................................................................................
## Purpose: To create a pretty figure for the manuscript
##
## Notes: 
## .................................................................................
source("R/00-functions_basic.R")
source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(raster)
library(srvyr) #wrap the survey package in dplyr syntax
library(RColorBrewer)


#............................................................
# read in data
#...........................................................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
dcdr <- readxl::read_excel(path = "model_datamaps/sub_DECODER_covariate_map_v3.xlsx", sheet = 1) %>% 
  dplyr::mutate(risk_factor_raw = ifelse(is.na(risk_factor_raw), "n", risk_factor_raw),
                risk_factor_model = ifelse(is.na(risk_factor_model), "n", risk_factor_model))
dtsrvy <- makecd2013survey(survey = dt)
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")
genosf <- ge %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2])
sf::st_geometry(genosf) <- NULL
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")

#............................................................
# read in results
#...........................................................
# point
point_prev <- readRDS("results/Pointest_map_prev_PlotObj.RDS")  + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))
# prov
prov_prev <- readRDS("results/ProvMap_postMeans_PlotObj.RDS") + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))
prov_se <- readRDS("results/ProvMap_postMeans_PlotObj.RDS")  + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))
prov_se <- readRDS("results/ProvMap_postSE_PlotObj.RDS")  + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))
# clust
clust_prev <- readRDS("results/ClustPrevMap_postMeans_PlotObj.RDS")  + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))
clust_se <- readRDS("results/ClustPrevMap_postSE_PlotObj.RDS")  + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))

#............................................................
# get binomial exact SE for clusters
#...........................................................
# cluster-level prevalence, because everyone is weighted the same in the cluster, don't use weights here
# note, the numerators will be slightly different (e.g. N) but the denomminators adjust for this
# Going to report whole numbers/unadjusted for clusters
clst <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(n = srvyr::survey_total(count),
                   pv18sn = srvyr::survey_total(pv18s)) %>% 
  dplyr::mutate(
    pv18sbinomtest = purrr::map2(pv18sn, n, function(x,y)
    {binom.test(x = round(x), n = round(y))}),
    pv18sprev = pv18sn/n,
    pvconfint = purrr::map(pv18sbinomtest, "conf.int"),
    pv18sL95 = purrr::map_dbl(pvconfint, 1),
    pv18sU95 = purrr::map_dbl(pvconfint, 2),
    se = (pv18sU95 - pv18sL95)/3.92)

point_se <- dplyr::left_join(clst, genosf, by = "hv001") %>% 
  ggplot() + 
  geom_sf(data = DRCprov, color = "#737373", fill = "#525252") +
  geom_point(aes(x = longnum, y = latnum, color = se)) +
  scale_color_distiller("Standard Error", type = "div", palette = "RdYlBu",
                       na.value = NA) + 
  theme_void() + 
  theme(legend.position = "bottom",
        legend.text = element_text(angle = 45, hjust = 0.5, vjust = 0.5, face = "bold"),
        plot.margin = margin(0, 0, 1, 0.1, "cm"))

jpeg("~/Desktop/temp2.jpg", width = 11, height = 8, res = 500, units = "in")
cowplot::plot_grid(point_prev, prov_prev, clust_prev,
                   labels = c("(A)", "(B)", "(C)"), nrow = 1)
graphics.off()

jpeg("~/Desktop/temp3.jpg", width = 11, height = 8, res = 500, units = "in")
cowplot::plot_grid(point_se, prov_se, clust_se,
                   labels = c("(A)", "(B)", "(C)"), nrow = 1,
                   rel_heights = c(0.1, 0.45, 0.45))
graphics.off()


#............................................................
# give me distribution of cluster prev 
#...........................................................
clstmodel <- readRDS("analyses/07-spatial_prediction/prevmap_predictions/covar_predictions.RDS")
clstmodel.fitted <- expit(clstmodel$samples)
clstmodel.coords <- clstmodel$grid.pred
clstmodel.fitted.postmeans <- apply(clstmodel.fitted, 2, mean)
summary(clstmodel.fitted.postmeans)
sum(clstmodel.fitted.postmeans < 0.0296 & clstmodel.fitted.postmeans > 0, na.rm  = T)
sum(!is.na(clstmodel.fitted.postmeans))
