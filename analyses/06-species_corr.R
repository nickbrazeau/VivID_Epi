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
  mutate(plsmdprevscale = scale(logitplsmdprev, center = T, scale = T))
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





#...........................................................................................
# Perform a Monte Carlo Sampling Scheme to Assess for Interference at the national level
#...........................................................................................

spclftvr <- data.frame(
  pfldh_fctb = c("falpos", "falpos", "falpos", "falpos",   "falneg", "falneg", "falneg"),
  pv18s_fctb = c("vivneg", "vivpos", "vivneg", "vivpos",   "vivpos", "vivpos", "vivneg"),
  po18s_fctb = c("ovneg",  "ovneg",  "ovpos",  "ovpos",    "ovpos",  "ovneg",  "ovpos"),
  species =    c("Pfmono", "Pf/Pv",  "Pf/Po",  "Pf/Pv/Po", "Pv/Po",  "Pvmono", "Pomono")
)


dt$malariastatus <- as.numeric( apply(dt[,c("pfldh", "pv18s", "po18s")], 1, sum) >= 1 ) # malaria cases

df <- dt %>% 
  dplyr::filter(malariastatus == 1 ) %>% 
  dplyr::left_join(x=., y = spclftvr) %>% 
  dplyr::mutate(species = factor(species))

test_prop <- c("Pf"=sum(df$pfldh),
               "Po"=sum(df$po18s),
               "Pv"=sum(df$pv18s))

#test_prop_prob <- test_prop/sum(test_prop)
test_prop_prob <- test_prop/nrow(df)

probs = c("Pfmono"   =  test_prop_prob["Pf"] * (1-test_prop_prob["Po"]) * (1 - test_prop_prob["Pv"]),
          "Pf/Pv"    =  test_prop_prob["Pf"] * (1-test_prop_prob["Po"]) * test_prop_prob["Pv"],  
          "Pf/Po"    =  test_prop_prob["Pf"] * test_prop_prob["Po"] * (1-test_prop_prob["Pv"]),  
          "Pf/Pv/Po" =  test_prop_prob["Pf"] * test_prop_prob["Po"] * test_prop_prob["Pv"], 
          "Pv/Po"    =  (1-test_prop_prob["Pf"]) * test_prop_prob["Po"] * test_prop_prob["Pv"],  
          "Pvmono"   =  (1-test_prop_prob["Pf"]) * (1-test_prop_prob["Po"]) * test_prop_prob["Pv"], 
          "Pomono"   =  (1-test_prop_prob["Pf"]) * test_prop_prob["Po"] * (1-test_prop_prob["Pv"])
)
names(probs) <- gsub(".Pf", "", names(probs)) # weird R behavior for vector names -- whoops

probs <- probs[sort(names(probs))] # to cooperate with table
chi <- chisq.test(table(df$species),
                  p = probs/sum(probs))
chi$expected
chi$observed
chi$residuals

siminfxn <- function(smplsz = 1e3, gamma, species){
  ret <- sample(x = species, size = smplsz, replace = T, prob = gamma)
  return(ret)
}

reps <- 1e4
exp <- replicate(reps, siminfxn(smplsz = nrow(df), gamma = probs, 
                                species = names(probs))) %>% 
  as.data.frame(.) %>% 
  magrittr::set_colnames(., paste0("sim", seq(1:ncol(.)))) %>% 
  tidyr::gather(., key = "simnum", value = "species") %>% 
  dplyr::mutate(species = factor(species))

# group_by does not respect the 0 levels for factors, can just use xtabs for simple cae
exp <- as.data.frame( xtabs(~species + simnum, data = out) )

obs <- as.data.frame( xtabs(~species, data = df) )
# plot it
ggplot() +
  geom_bar(data = exp, aes(x = Freq, y=..count..)) +
  geom_vline(data = obs, aes(xintercept = Freq), color = "red") +
  facet_wrap(~species, scales="free")



#...........................................................................................
# At the Household level are Pf and Pv interferring
#...........................................................................................
# dtsrvy %>% 
#   dplyr::group_by(houseid) %>% 
#   dplyr::mutate(count = 1) %>% 
#   dplyr::summarise(n = srvyr::survey_total(count), 
#                    pfldhn = srvyr::survey_total(pfldh, na.rm = T), 
#                    pfldhprev = srvyr::survey_mean(pfldh, na.rm = T, vartype = c("se", "ci"), level = 0.95),
#                    pv18sn = srvyr::survey_total(pv18s, na.rm = T), 
#                    pv18sprev = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95)
#                    ) %>%
#   ggplot(.) +
#   geom_point(aes(x=pfldhprev, y = pv18sprev))

pvptonpfrst <- dtsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    n = srvyr::survey_total(count),
    pvcases = srvyr::survey_total(pv18s),
    pvprev = srvyr::survey_mean(pv18s)) %>% 
  dplyr::left_join(x=., y = ge)

pos <- pvptonpfrst %>% 
  dplyr::filter(pvprev > 0)
neg <- pvptonpfrst %>% 
  dplyr::filter(pvprev == 0)


prevmaprasterplots[[1]] + list(
  geom_sf(data = DRCprov, fill = "NA"),
  geom_jitter(data = neg, aes(x=longnum, y=latnum, size = n), shape = 4, show.legend = F, colour = "#000000", alpha = 0.8),
  geom_jitter(data = pos, aes(x=longnum, y=latnum, colour = pvprev, size = n), alpha = 0.8),
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26") 
)




#...........................................................................................
# Proportion of Mono-infections
#...........................................................................................
pvmono <- dtsrvy %>% 
  dplyr::filter(pv18s == 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
                   pvcases = srvyr::survey_total(pv18s, na.rm = T),
                   pvmono = srvyr::survey_total(pv18s == 1 & pfldh ==0, na.rm = T)) 
                     
pvmono$pvmonoprop <- pvmono$pvmono/pvmono$pvcases
pvmono <- pvmono %>% 
  left_join(., y=ge) 

pvmonomap <- pvmono %>% 
  ggplot(.) +
  geom_sf(data = DRCprov, fill = "NA") +
  geom_point(aes(x=longnum, y=latnum, color = pvmonoprop)) +
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26") 

plot(pvmonomap)

prevmaprasterplots[[1]] + list(
    geom_sf(data = DRCprov, fill = "NA"),
    geom_point(data = pvmono, aes(x=longnum, y=latnum, color = pvmonoprop)),
    scale_colour_gradient(low = "#2b8cbe", high = "#de2d26") 
    )



#...........................................................................................
# Remove singleton clusters
#...........................................................................................
pvnosing <- dtsrvy %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    pvcases = srvyr::survey_total(pv18s),
    pvprev = srvyr::survey_mean(pv18s)) %>% 
  dplyr::filter(pvcases > 2) %>% 
  left_join(., y=ge) %>% 
  ggplot(.) +
  geom_sf(data = DRCprov, fill = "NA") +
  geom_point(aes(x=longnum, y=latnum, color = pvprev)) +
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26") 
plot(pvnosing)


pvnosing <- dt %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
    pvcases = sum(pv18s),
    pvprev = mean(pv18s)) %>% 
  dplyr::filter(pvcases > 1) %>% 
  left_join(., y=ge) %>% 
  ggplot(.) +
  geom_sf(data = DRCprov, fill = "NA") +
  geom_point(aes(x=longnum, y=latnum, color = pvprev)) +
  scale_colour_gradient(low = "#2b8cbe", high = "#de2d26") 
plot(pvnosing)



