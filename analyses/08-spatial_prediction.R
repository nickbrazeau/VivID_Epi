#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create a spatial prediction model
#----------------------------------------------------------------------------------------------------
source("~/Documents/GitHub/VivID_Epi/R/00-functions_maps.R")

#......................
# Import Data
#......................
dt <- readRDS("~/Documents/GitHub/VivID_Epi/data/derived_data/vividepi_recode.rds")
options(survey.lonely.psu="certainty")
dtsrvy <- dt %>% srvyr::as_survey_design(ids = hv001, strata = hv023, weights = hv005_wi)

#spatial from the DHS -- these are cluster level vars
ge <- sf::st_as_sf(readRDS(file = "~/Documents/GitHub/VivID_Epi/data/raw_data/dhsdata/datasets/CDGE61FL.rds"))
colnames(ge) <- tolower(colnames(ge))
ge$adm1name <- gsub(" ", "-", ge$adm1name) # for some reason some of the char (like Kongo Central, lack a -), need this to match GADM
ge$adm1name <- gsub("Tanganyka", "Tanganyika", ge$adm1name) # DHS misspelled this province
# remove clusters that were missing from the DHS, see readme
ge <- ge %>% 
  dplyr::rename(hv001 = dhsclust) # for easier merge with PR


#......................
# Summarize by Cluster
#......................
pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001") %>% 
  dplyr::left_join(x=., y = ge)
pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001") %>% 
  dplyr::left_join(x=., y = ge)
po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001") %>% 
  dplyr::left_join(x=., y = ge)



#----------------------------------------------------------------------------------------------------
# Smoothed Guassian Maps
#----------------------------------------------------------------------------------------------------

# bind those to a tibble
pr <- dplyr::bind_rows(pfldhclust, pv18sclust, po18sclust) %>% 
  dplyr::group_by(plsmdmspec) %>% 
  tidyr::nest()


#.............................
# get prev rasters for individual levels
#..............................
# polybb <- osmdata::getbb("Democratic Republic of the Congo", featuretype = "country",  format_out = 'polygon')
poly <- cbind(c(17,32,32,12,12), c(-14,-14,6,6,-14)) 
grid.pred <- splancs::gridpts(poly, xs=0.1, ys=0.1)
colnames(grid.pred) <- c("long","lat")

pr$prevrasters <- map(pr$data, 
                      fit_pred_spMLE, outcome = "logitplsmdprev", covar = "1", 
                      long_var = "longnum", lat_var = "latnum",
                      grid.pred = grid.pred, kappa = 0.5, 
                      pred.reps = 10)

pr$prevrasterspred <- purrr::map(pr$prevrasters, "pred")



#.............................
# plot prev rasters
#..............................
prevmaprasterplots <- lapply(pr$prevrasterspred,
                             prevmaprasterplotter, smoothfct = rep(7,3))
prevmaprasterplots <- map(prevmaprasterplots, function(x){return(x + prettybasemap_nodrc)})


jpeg(file = "~/Documents/GitHub/VivID_Epi/figures/04-guassian-prev-maps.jpg", width = 11, height = 8, units="in", res=300)
gridExtra::grid.arrange(
  prevmaprasterplots[[1]], 
  prevmaprasterplots[[2]],
  prevmaprasterplots[[3]],
  nrow=1, 
  top=grid::textGrob("Smoothed Prevalence by Species in CD2013 DHS", gp=grid::gpar(fontsize=15, fontfamily = "Arial", fontface = "bold"))) 
graphics.off()



#.............................
# prev clusters
#..............................

#..........................
# spatial kernel densities
#..........................
cases <- dt[dt$pv18s == 1, ]

casedens <- MASS::kde2d(x = cases$longnum, 
                 y = cases$latnum,
                 n=1e3,
            lims = c(bb[1,], bb[2,]))

contour(casedens)


controls <- dt[dt$pv18s == 0, ]
condens <- MASS::kde2d(x = controls$longnum, 
                    y = controls$latnum,
                    n=1e3,
                    lims = c(bb[1,], bb[2,]))

contour(condens)



