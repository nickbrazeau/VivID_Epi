#----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate for signals that shouldn't be there
# This is to better understand spatial randomnes of Pv
#----------------------------------------------------------------------------------------------------
library(tidyverse)
library(RColorBrewer)
library(raster)
library(gstat)
load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
load("~/Documents/GitHub/VivID_Epi/data/vividmaps_small.rda")
load("~/Documents/GitHub/VivID_Epi/data/vividmaps_large.rda")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_basic.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_glms.R")
source("~/Documents/GitHub/VivID_Epi/R/00-functions_guassmap.R")
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum", "geometry")]
#----------------------------------------------------------------------------------------------------
# Check for Contamination
#----------------------------------------------------------------------------------------------------
#........................
# Original Plate nums
#........................
plts <- dt %>% 
  group_by(original_platemnum, hv001) %>% 
  summarise(longnum = mean(longnum), latnum = mean(latnum))

pltnm <- dt %>% 
  group_by(original_platemnum) %>% 
  summarise(n = n(), 
            pvcount = sum(pv18s, na.rm=T),
            pvmean = mean(pv18s, na.rm=T),
            pvctmean = mean(pv18sct_cont, na.rm = T)) %>% 
  filter(pvmean > 0) 
  

pltnm %>% 
  DT::datatable(., extensions='Buttons',
                options = list(
                  searching = T,
                  pageLength = 20,
                  dom = 'Bfrtip', 
                  buttons = c('csv')))
# https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r              

n <- length(levels(factor(plts$hv001)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- c(col_vector, rainbow(7))

ggplot() +
  geom_sf(data = DRCprov, fill = "#d9d9d9") +
  geom_point(data = plts, aes(x=longnum, y=latnum, colour = factor(original_platemnum))) +
#  scale_color_manual(values = col_vector) +
  theme(legend.position = "none")


ggplot() +
  geom_sf(data = DRCprov, fill = "#d9d9d9") +
  ggrepel::geom_text_repel(data = clustgeom, aes(x=longnum, y=latnum, label = hv001)) +
  theme(legend.position = "none")


hvmaps <- dt %>% 
  group_by(original_platemnum) %>% 
  mutate(hvcat = paste(levels(factor(hv001)), collapse = "-")) %>% 
  select(c("original_platemnum", "hvcat")) 
hvmaps <- hvmaps[!duplicated(hvmaps),]

p1 <- dt %>% 
  group_by(original_platemnum) %>% 
  summarise(n = n(), 
            pvcount = sum(pv18s, na.rm=T),
            pvmean = mean(pv18s, na.rm=T),
            pvctmean = mean(pv18sct_cont, na.rm = T)) %>% 
  filter(pvmean > 0) %>% 
  left_join(x=., y=hvmaps, by = "original_platemnum") %>% 
  mutate(original_platemnum_num = as.numeric(factor(original_platemnum))) %>% 
  ggplot() +
  geom_point(aes(x=original_platemnum_num, y=pvmean, pltnm = original_platemnum,
                 hvcat = hvcat)) 
plotly::ggplotly(p1)

#........................
# Are plates overly contaminated?
#........................
platmaps <- dt %>% 
  group_by(hv001) %>% 
  mutate(original_platemnum_cat = paste(levels(factor(original_platemnum)), collapse = "-")) %>% 
  select(c("hv001", "original_platemnum_cat")) 
platmaps <- platmaps[!duplicated(platmaps),]

p2 <- dt %>% 
  group_by(hv001) %>%   
   summarise(n = n(), 
            pvcount = sum(pv18s, na.rm=T),
            pvmean = mean(pv18s, na.rm=T),
            pvctmean = mean(pv18sct_cont, na.rm = T)) %>% 
  filter(pvmean > 0) %>% 
  left_join(x=., y=platmaps, by = "hv001") %>% 
  ggplot() +
  geom_point(aes(x=hv001, y=pvmean, platnum = original_platemnum_cat)) 
plotly::ggplotly(p2)

pltnm %>% 
  filter(original_platemnum == "M061") %>% 
  ggplot() +
  geom_sf(data = DRCprov, fill = "#d9d9d9") +
  ggrepel::geom_text_repel(aes(x=longnum, y=latnum, label = hv001)) +
  theme(legend.position = "none")



#----------------------------------------------------------------------------------------------------
# Sampling Dates Distributions
#----------------------------------------------------------------------------------------------------
datesamp <- dt %>% 
  dplyr::group_by(hv001) %>%
  dplyr::mutate(time = ifelse(length(unique(hvyrmnth_fctm)) == 1, 
                              paste(unique(hvyrmnth_fctm)),
                              paste0(unique(hvyrmnth_fctm), collapse="/"))
  ) %>% 
  dplyr::summarise(n = n(), 
                   timen = length(unique(hvyrmnth_fctm)),
                   latnum = mean(latnum), 
                   longnum = mean(longnum),
                   time = unique(time)
  ) 
datesampmap <- datesamp %>% 
  ggplot() +
  geom_sf(data = DRCprov, fill = "#d9d9d9") +
  geom_point(data = datesamp, aes(x=longnum, y=latnum, color = factor(time)),
             size = 1.2, alpha = 0.8) +
  theme(legend.position = "bottom",
        plot.background = element_blank())



jpeg("~/Documents/GitHub/VivID_Epi/figures/07-clstr_collection_times.jpg", width = 8, height = 8, res = 400, units = "in")
plot(datesampmap)
graphics.off()





#----------------------------------------------------------------------------------------------------
# Quick Sensitivity Analysis 
#----------------------------------------------------------------------------------------------------
#........................
# recode and rerun
#........................
# Have EMPIRICAL reason to believe that this 2 CT off...
dt$pv18s_sens <- ifelse(dt$pv18sct_cont < 42 & !is.na(dt$pv18sct_cont), 1, 0)
pv18sprov_sens <- prev_point_est_summarizer(data = dt, maplvl = adm1name, plsmdmspec = pv18s_sens)  %>% 
  dplyr::rename(plsmdprev_sens = plsmdprev)
pv18sclust_sens <- prev_point_est_summarizer(data = dt, maplvl = hv001, plsmdmspec = pv18s_sens) %>% 
  magrittr::set_colnames( c(colnames(.)[colnames(.) %in% c("hv001", "n")],
                            paste0(colnames(.)[!colnames(.) %in% c("hv001", "n")], "_sens" ) ) )
# get unweighted
unweighted_pv_sen <- dt %>% 
  group_by(hv001) %>% 
  summarise(unweighted_count_sens = sum(pv18s_sens, na.rm = T),
            unweighted_prev_sens = mean(pv18s_sens, na.rm = T))

pv18sclust_sens <- left_join(pv18sclust_sens, unweighted_pv_sen)

#......................
# Plot Cases
#......................
casemapplotter <- function(data, plsmdmspec, prev_var){
  
  prev_var <- dplyr::enquo(prev_var)
  
  # Set some colors ; took this from here https://rjbioinformatics.com/2016/07/10/creating-color-palettes-in-r/ ; Here is a fancy color palette inspired by http://www.colbyimaging.com/wiki/statistics/color-bars
  clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum")]
  data <- inner_join(data, clustgeom, by = "hv001")
  pos <- data %>% 
    dplyr::filter(!!prev_var > 0)
  neg <- data %>% 
    dplyr::filter(!!prev_var == 0)
  
  ret <- ggplot() + 
    geom_sf(data = DRCprov) +
    geom_jitter(data = neg, aes(x=longnum, y=latnum, size = n), shape = 4, show.legend = F, colour = "#377eb8") +
    geom_point(data = pos, aes(x=longnum, y=latnum, colour = plsmdprev_sens, size = n), alpha = 0.4) +
    scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
    scale_size(guide = 'none') +
    ggtitle(paste(plsmdmspec)) +
    coord_sf(datum=NA) + # to get rid of gridlines
    vivid_theme +
    theme(axis.text = element_blank(),
          axis.line = element_blank(), 
          axis.title = element_blank(),
          legend.position = "bottom")
  
  return(ret)
  
}

casemapplotter(pv18sclust_sens, "pv18s", prev_var = plsmdprev_sens)
hist(pv18sclust_sens$plsmdprev[pv18sclust_sens$plsmdprev_sens != 0])
min(pv18sclust_sens$plsmdprev[pv18sclust_sens$plsmdprev_sens != 0])
sum(pv18sclust_sens$plsmdprev_sens > 0)


# ORIGINAL 

orig <- prev_point_est_summarizer(data = dt, maplvl = hv001, plsmdmspec = pv18s) %>% 
  magrittr::set_colnames( c(colnames(.)[colnames(.) %in% c("hv001", "n")],
                    paste0(colnames(.)[!colnames(.) %in% c("hv001", "n")], "_orig" ) ) )
casemapplotter(orig, "pv18s", prev_var = plsmdprev_orig)

# ORIGINAL unweighted
unweighted_pv_orig <- dt %>% 
  group_by(hv001) %>% 
  summarise(unweighted_count_orig = sum(pv18s, na.rm = T),
            unweighted_prev_orig = mean(pv18s, na.rm = T))

orig <- left_join(orig, unweighted_pv_orig)


hist(orig$plsmdprev_orig[orig$plsmdprev_orig != 0])
min(orig$plsmdprev_orig[orig$plsmdprev_orig != 0])
sum(orig$plsmdprev_orig > 0)

jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/07-vivax-prev_sens.jpg", width = 11, height = 8, units="in", res=300)
casemapplotter(pv18sclust_sens, "pv18s_sense")
graphics.off()

# which clusters are lost
anti <- dplyr::anti_join(x=orig[orig$plsmdprev != 0, ], 
                         y=pv18sclust_sens[pv18sclust_sens$plsmdprev != 0, ], 
                         by = "hv001")


corr <- inner_join(pv18sclust_sens, orig, by = c("hv001", "n"))
nrow(corr) == nrow(orig)

ggplot() +
  geom_point(data = corr, aes(x=plsmdprev_sens, y=plsmdprev_orig, size = n), alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color="red", alpha = 0.5) + 
  ggtitle("Cluster Prevalences") + xlab("CT < 42") + ylab("CT < 45 & Snounou") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
write.csv(corr, file = "reports/report_obj/correlation_of_Pv_from_Sensitivity_Analysis.csv",
          row.names = F, quote = F)

corr %>% 
  filter(unweighted_count_orig > 0) %>% 
  ggplot() + 
  geom_bar(data = , aes(x=unweighted_count_orig)) +
  scale_x_continuous("Unweighted Count of Cases within a Cluster", breaks = seq(1:10)) +
  ylab("Count of the Count (i.e. Count of Number of Cases)") + 
  ggtitle("Distribution of the Number of Cases in a Cluster") +
  theme(plot.title = element_text(hjust=0.5))

#.................
# Parametric, Bivariate Analysis
# Odds Ratios with Pv as the outcome
#.................
# see 03-Uni_Bivar_Analyses_SENS



#----------------------------------------------------------------------------------------------------
# Looking at Singleton Clusters
#----------------------------------------------------------------------------------------------------

clstbarcodect <- dt[, c("hv001", "pv18sct_cont", "hivrecode_barcode")]

clstbarcodectcorr <- left_join(clstbarcodect, corr, by = "hv001") %>% 
  filter(!is.na(pv18sct_cont)) %>% 
  filter(unweighted_count_sens != 0 &  unweighted_count_orig != 0) # didn't pass snounou


clstbarcodectcorr %>% 
  mutate_if(is.numeric, round, 2) %>% 
  DT::datatable(., extensions='Buttons',
                options = list(
                  searching = T,
                  pageLength = 20,
                  dom = 'Bfrtip', 
                  buttons = c('csv')))

# are cts different?
clstbarcodectcorr %>% 
  ggplot(data = .) +
  ggridges::geom_density_ridges(aes(x=pv18sct_cont, y = factor(unweighted_count_orig))) +
  ylab("Count of Cases per cluster")


# remove singletons and map 
corr %>% 
  dplyr::mutate(plsmdprev_orig_nosing = ifelse(unweighted_count_orig == 1, 0, plsmdprev_orig)) %>% 
  casemapplotter(., "pv18s", prev_var = plsmdprev_orig_nosing)

# remove singletons/doubles and map 
corr %>% 
  dplyr::mutate(plsmdprev_orig_nosing = ifelse(unweighted_count_orig %in% c(1,2), 0, plsmdprev_orig)) %>% 
  casemapplotter(., "pv18s", prev_var = plsmdprev_orig_nosing)

# remove singletons/doubles/tripletons and map 
corr %>% 
  dplyr::mutate(plsmdprev_orig_nosing = ifelse(unweighted_count_orig %in% c(1,2,3), 0, plsmdprev_orig)) %>% 
  casemapplotter(., "pv18s", prev_var = plsmdprev_orig_nosing)

# repel for positive clusters 
pt <- corr %>% 
  casemapplotter(., "pv18s", prev_var = plsmdprev_orig) 

repel <- corr %>% 
  dplyr::filter(unweighted_count_orig >= 1) %>% 
  left_join(x=., y=clustgeom) 

pt +
  ggrepel::geom_label_repel(data = repel, aes(x=longnum, y = latnum, label = unweighted_count_orig))

#----------------------------------------------------------------------------------------------------
# QUICK Interpolation for Temperature and Precip
#----------------------------------------------------------------------------------------------------
# https://mgimond.github.io/Spatial/interpolation-in-r.html
tempprecip <- dt %>% 
  group_by(hv001) %>% 
  summarise(clsttemp = mean(mean_temperature_2015_cont),
            clstsdtemp = sd(mean_temperature_2015_cont),
            clstrain = mean(rainfall_2015_cont),
            clstsdrain = sd(rainfall_2015_cont)) %>% 
  left_join(., y=clustgeom, by = "hv001")

sum(apply(tempprecip[,c(3,5)], 1, function(x){return(x > 0)}), na.rm = T)


# imputing for mean
tempprecip$clsttemp[is.na(tempprecip$clsttemp)] <- mean(tempprecip$clsttemp, na.rm=T)

#---------------------------------------
# Temperature
#---------------------------------------

#.....................
# make raster and boundaries
#....................
polybb <- getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long", "lat")
pos <- sf::as_Spatial(sf::st_as_sf(polybb))

st_sf(bb = 1:2, geom = st_sfc( st_point(bb[,1]), st_point(bb[,2]) ) )


grd              <- as.data.frame(spsample(grid.pred, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(pos)

temp.idw <- gstat::idw(tempprecip$clsttemp ~ 1, pos, 
                    newdata = grd, idp = 2.0)
# Convert to raster object then clip to DRC
temp.r       <- raster(temp.idw)
temp.r.m     <- mask(temp.r, as_Spatial(DRCprov))

# First smooth and then make the temperature raster plot
temp.r.m <- focal(temp.r.m, w=matrix(1,
                                       nrow=5,
                                       ncol=5), mean)
temp.r.m.pts  <-  raster::rasterToPoints(temp.r.m)
temp.r.m.pts.df <-  data.frame(lon = temp.r.m.pts[,1], 
                               lat = temp.r.m.pts[,2], 
                               temp = temp.r.m.pts[,3])
temp.r.m.pts.m.plot <- ggplot() + 
  geom_raster(data = temp.r.m.pts.df, aes(lon, lat, fill = temp), alpha = 0.8) +
  scale_fill_gradient2("Temperature", low = "#313695", mid = "#ffffbf", high = "#a50026",
                       midpoint = 21) +
  prettybasemap_nodrc 



jpeg("~/Documents/GitHub/VivID_Epi/figures/07-temperature-idw.jpg", 
     height = 8, width=8, units = "in", res=400)
plot(temp.r.m.pts.m.plot)
graphics.off()


#---------------------------------------
# Precipation
#---------------------------------------
precip.idw <- gstat::idw(tempprecip$clstrain ~ 1, pos, 
                       newdata = grd, idp = 2.0)
# Convert to raster object then clip to DRC
precip.r       <- raster(precip.idw)
precip.r.m     <- mask(precip.r, as_Spatial(DRCprov))

# First smooth and then make the temperature raster plot
precip.r.m <- focal(precip.r.m, w=matrix(1,
                                     nrow=5,
                                     ncol=5), mean)
precip.r.m.pts  <-  raster::rasterToPoints(precip.r.m)
precip.r.m.pts.df <-  data.frame(lon = precip.r.m.pts[,1], 
                               lat = precip.r.m.pts[,2], 
                               precip = precip.r.m.pts[,3])
precip.r.m.pts.m.plot <- ggplot() + 
  geom_raster(data = precip.r.m.pts.df, aes(lon, lat, fill = precip), alpha = 0.8) +
  scale_fill_gradient2("Precipitation", low = "#313695", mid = "#ffffbf", high = "#a50026",
                       midpoint = 1520)+
  prettybasemap_nodrc 



jpeg("~/Documents/GitHub/VivID_Epi/figures/07-precip-idw.jpg", 
     height = 8, width=8, units = "in", res=400)
plot(precip.r.m.pts.m.plot)
graphics.off()


