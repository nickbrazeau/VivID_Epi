#----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate if the pv and pf data are correlated
# considering both space and other epi variables
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(srvyr)
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
load("~/Documents/GitHub/VivID_Epi/data/04-basic_mapping_data.rda")
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "hv025_fctb", "latnum", "longnum", "geometry")]

#...................................
# Loess Plots based on Summary 
#...................................
pfclst <- mp$data[[4]]
colnames(pfclst)[6] <- "Pf_prev"

pvclst <- mp$data[[5]]
colnames(pvclst)[6] <- "Pv_prev"


loesspvpf <- dplyr::bind_cols(pfclst, pvclst) %>% 
  ggplot(aes(x=Pf_prev, y=Pv_prev, size = n)) +
  geom_point(aes(colour = hv001), show.legend = F) +
  geom_smooth(aes(x=Pf_prev, y=Pv_prev, weight = n), method="loess", se=F, colour = "red", show.legend = F) +
  ggtitle("Pfalciparum versus Pvivax Cluster Prevalence") + ylab("Pv Prevalence") + xlab("Pf Prevalence") + labs(caption = "Removed 0,0 Clusters") + 
  vivid_theme
plotly::ggplotly(loesspvpf)
jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/06-pv_vs_pf.jpg", width = 11, height = 8, units = "in", res=300)
plot(loesspvpf)
graphics.off()

loesspvpf_rrlurbn <- dplyr::bind_cols(pfclst, pvclst) %>% 
  left_join(., y=clustgeom) %>% 
  ggplot(aes(x=Pf_prev, y=Pv_prev, size = n, group = hv025_fctb)) +
  geom_point(aes(colour = hv001), show.legend = F) +
  geom_smooth(aes(x=Pf_prev, y=Pv_prev, weight = n), method="loess", se=F, colour = "red", show.legend = F) +
  ggtitle("Pfalciparum versus Pvivax Cluster Prevalence") + ylab("Pv Prevalence") + xlab("Pf Prevalence") + labs(caption = "Removed 0,0 Clusters") + 
  facet_wrap(~hv025_fctb) +
  vivid_theme
plotly::ggplotly(loesspvpf_rrlurbn)

jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/06-pv_vs_pf_byurban.jpg", width = 11, height = 8, units = "in", res=300)
plot(loesspvpf_rrlurbn)
graphics.off()





#...................................
# Distance Matrix
#...................................
pfclst <- mp$data[[4]]
pvclst <- mp$data[[5]]





#...................................
# Comparison of Pf and Pv Age dist
#...................................
pvage <- dt %>% 
  dplyr::mutate(count = 1) %>% 
  srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
  dplyr::group_by(hv105) %>% 
  dplyr::summarise(plsmdn = srvyr::survey_total(count), 
                   plsmd = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95))
pvage$species <- "vivax"

pfage <- dt %>% 
  dplyr::mutate(count = 1) %>% 
  srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
  dplyr::group_by(hv105) %>% 
  dplyr::summarise(plsmdn = srvyr::survey_total(count), 
                   plsmd = srvyr::survey_mean(pfldh, na.rm = T, vartype = c("se", "ci"), level = 0.95))

pfage$species <- "falciparum"

plotObj <- dplyr::bind_rows(pvage, pfage) %>% 
  ggplot() + 
  geom_line(aes(x=hv105, y=plsmd, color = species)) + 
  geom_ribbon(aes(x=hv105, ymin=plsmd_low, ymax=plsmd_upp, fill = species), alpha=0.5) +
  geom_point(aes(x=hv105, y=plsmd, size=plsmdn, color = species), alpha=0.5, show.legend=F) +
  scale_fill_manual("Malara Species", values = c("#a50026", "#313695"), labels = c("Pfalciparum", "Pvivax")) +
  scale_colour_manual(values = c("#a50026", "#313695"), guide =F) +
  ggtitle("Malaria Prevalence by Age") +
  xlab("Age") + ylab("Malaria Prevalence") + 
  vivid_theme



