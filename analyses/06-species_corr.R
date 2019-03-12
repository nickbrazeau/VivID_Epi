#----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate if the pv and pf data are correlated
# considering both space and other epi variables
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(srvyr)
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R") 
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R") 

#......................
# Import Data
#......................
# epi data
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
options(survey.lonely.psu="certainty")
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, strata = hv023, weights = hv005_wi)

# base maps
load("~/Documents/GitHub/VivID_Epi/data/map_bases/vivid_maps_bases.rda")

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


#......................
# Cluster and Prov Level Data
#......................
pfldhprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pfldh) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "adm1name")
pv18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pv18s)  %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "adm1name")
po18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = po18s) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "adm1name")

pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001")
pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001")
po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001")



# bind those to a tibble
mp <- dplyr::bind_rows(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust) %>% 
  dplyr::group_by(plsmdmspec, maplvl) %>% 
  tidyr::nest()


# this awful hack becuase of this issue https://github.com/tidyverse/dplyr/issues/3483
# we are going down the rabbit hole just to try and make this stupid survey and purr package work. fine for now but return
mp$data <- lapply(list(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust), function(x) return(x))


#...................................
# Are CT values correlated?
#...................................
summary( lm(dt$pv18sct_cont_ind_log_scale ~ dt$pfldhct_cont_ind_log_scale) )
plot(dt$pv18sct_cont_ind_log_scale ~ dt$pfldhct_cont_ind_log_scale)

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

#jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/06-pv_vs_pf.jpg", width = 11, height = 8, units = "in", res=300)
#plot(loesspvpf)
#graphics.off()

loesspvpf_rrlurbn <- dplyr::bind_cols(pfclst, pvclst) %>% 
  left_join(., y=ge) %>% 
  ggplot(aes(x=Pf_prev, y=Pv_prev, size = n, group = urban_rura)) +
  geom_point(aes(colour = hv001), show.legend = F) +
  geom_smooth(aes(x=Pf_prev, y=Pv_prev, weight = n), method="loess", se=F, colour = "red", show.legend = F) +
  ggtitle("Pfalciparum versus Pvivax Cluster Prevalence") + ylab("Pv Prevalence") + xlab("Pf Prevalence") + labs(caption = "Removed 0,0 Clusters") + 
  facet_wrap(~urban_rura) +
  vivid_theme
plotly::ggplotly(loesspvpf_rrlurbn)

# jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/06-pv_vs_pf_byurban.jpg", width = 11, height = 8, units = "in", res=300)
# plot(loesspvpf_rrlurbn)
# graphics.off()


#...................................
# Kids RDT and Micro
#...................................

loessrdtpf_rrlurbn <- left_join(pvclst, rdtmicro, by= "hv001") %>%
  dplyr::left_join(x=., y=ge) %>% 
  ggplot(aes(x=RDTprev, y=Pv_prev, size = n, group = urban_rura)) +
  geom_point(aes(colour = hv001), show.legend = F) +
  geom_smooth(aes(x=RDTprev, y=Pv_prev, weight = n), method="loess", se=F, colour = "red", show.legend = F) +
  ggtitle("Pfalciparum RDT versus Pvivax Cluster Prevalence") + ylab("Pv Prevalence") + xlab("Pf Prevalence") + labs(caption = "Removed 0,0 Clusters") + 
  facet_wrap(~urban_rura) +
  vivid_theme
plotly::ggplotly(loessrdtpf_rrlurbn)


loessmicropf_rrlurbn <- left_join(pvclst, rdtmicro, by= "hv001") %>%
  dplyr::left_join(x=., y=ge) %>% 
  ggplot(aes(x=microprev, y=Pv_prev, size = n, group = urban_rura)) +
  geom_point(aes(colour = hv001), show.legend = F) +
  geom_smooth(aes(x=RDTprev, y=Pv_prev, weight = n), method="loess", se=F, colour = "red", show.legend = F) +
  ggtitle("Pfalciparum RDT versus Pvivax Cluster Prevalence") + ylab("Pv Prevalence") + xlab("Pf Prevalence") + labs(caption = "Removed 0,0 Clusters") + 
  facet_wrap(~urban_rura) +
  vivid_theme
plotly::ggplotly(loessmicropf_rrlurbn)




#...................................
# Distance Matrix
#...................................
pfclst <- mp$data[[4]] %>% 
  mutate(plsmdprevscale = scale(logitplsmdprev, center = F, scale = T))
colnames(pfclst) <- gsub("plsmd", "pf", colnames(pfclst))
boxplot(pfclst$logitpfprev)

pvclst <- mp$data[[5]] %>% 
  mutate(plsmdprevscale = scale(logitplsmdprev, center = F, scale = T))
colnames(pvclst) <- gsub("plsmd", "pv", colnames(pvclst))
boxplot(pvclst$logitpvprev)


left_join(pfclst, pvclst) %>% 
  left_join(x=., y=ge) %>% 
  mutate(plsmdiff = pfprevscale - pvprevscale) %>% 
  ggplot(.) +
  geom_sf(data = DRCprov, fill = "NA") +
  geom_point(aes(x=longnum, y=latnum, color = plsmdiff)) +
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26")


pfpv %>% 
  dplyr::filter(pvprev == 0) %>% 
  ggplot(.) +
  geom_sf(data = DRCprov, fill = "NA") +
  geom_point(aes(x=longnum, y=latnum, color = pfprev)) +
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26")


clst0s <- pfpv %>% 
  dplyr::filter(pvprev == 0 & pfprev == 0) %>% 
  dplyr::select(hv001)

clst0s <- dt[dt$hv001 %in% clst0s$hv001, ]





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






pr$prevrasters <- map(pr$data, 
                      fit_pred_spMLE, outcome = "logitplsmdprev", covar = "1", 
                      long_var = "longnum", lat_var = "latnum",
                      grid.pred = grid.pred, kappa = 0.5, 
                      pred.reps = 1)

pr$prevrasterspred <- purrr::map(pr$prevrasters, "pred")




