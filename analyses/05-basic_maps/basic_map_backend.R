## .................................................................................
## Purpose: Calculate Moran's I once instead of every knit
##
## Notes: 
## .................................................................................
source("R/00-functions_basic.R")
source("R/00-functions_maps.R") 
library(tidyverse)
library(sf)
library(geosphere)
library(srvyr) #wrap the survey package in dplyr syntax
set.seed(48)
#............................................................
# Import Data
#...........................................................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")

# need to add back in geometry to dt
gegeoms <- ge %>% 
  dplyr::select(c("hv001", "geometry"))
dt <- dplyr::left_join(dt, gegeoms, by = "hv001")

dtsrvy <- makecd2013survey(survey = dt)
#............................................................
# Plasmodium Point Prevalence Maps (Province & Cluster)
#...........................................................
pfldhprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pfldh, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "adm1name")
pv18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = pv18s, adm1shp = DRCprov)  %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "adm1name")
po18sprov <- prev_point_est_summarizer(design = dtsrvy, maplvl = adm1name, plsmdmspec = po18s, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "adm1name")

pfldhclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pfldh, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pfldh", maplvl = "hv001")
pv18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = pv18s, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "pv18s", maplvl = "hv001")
po18sclust <- prev_point_est_summarizer(design = dtsrvy, maplvl = hv001, plsmdmspec = po18s, adm1shp = DRCprov) %>% 
  dplyr::mutate(plsmdmspec = "po18s", maplvl = "hv001")



# bind those to a tibble
mp <- dplyr::bind_rows(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust) %>% 
  dplyr::group_by(plsmdmspec, maplvl) %>% 
  tidyr::nest()


# this awful hack becuase of this issue https://github.com/tidyverse/dplyr/issues/3483
# we are going down the rabbit hole just to try and make this stupid survey and purr package work. fine for now but return
mp$data <- lapply(list(pfldhprov, pv18sprov, po18sprov, pfldhclust, pv18sclust, po18sclust), function(x) return(x))


#............................................................
# Moran's I for Prov
#...........................................................
#......................
# Make Adjacency Matrix for Pvf
#......................
drcprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
W.nb <- spdep::poly2nb(sf::as_Spatial(drcprov), row.names = drcprov$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)


moran_mc_wrapper <- function(data, listw){
  prev <- data$plsmdprev
  ret <- spdep::moran.mc(x = prev,
                         listw = spdep::mat2listw(listw),
                         alternative = "greater",
                         nsim = 1e5)
  return(ret)
}

#........................
# Moran's I for Prov
#........................
mp.prov <- mp %>% 
  dplyr::filter(maplvl == "adm1name") %>% 
  dplyr::ungroup()

mp.prov$moranIprov <- purrr::map(mp.prov$data, moran_mc_wrapper, listw = W)

moranIprov <- mp.prov %>% 
  dplyr::select(c("plsmdmspec", "moranIprov")) 

moranIprov$Istatistic <- purrr::map(moranIprov$moranIprov, "statistic")
moranIprov$pvalue <- purrr::map(moranIprov$moranIprov, "p.value")
#......................
# save out
#......................
saveRDS(moranIprov, "results/MoranI_prov.RDS")


#............................................................
# Moran's I for Cluster
#...........................................................
mp.clst <- mp %>% 
  dplyr::filter(maplvl == "hv001") %>% 
  dplyr::ungroup()
#........................
# Moran's I for Cluster Level
#........................
ge.nosf <- ge
sf::st_geometry(ge.nosf) <- NULL

# note, need to remove three clusters that were not screened for malaria
# one cluster due to no HIV testing in it(? DHS internal) and 
# two clusters lost due to contamination
maldhsclustlist <- unique( dt$hv001 )
ge.clst <- ge %>% 
  dplyr::filter(hv001 %in% maldhsclustlist)

# get GC Weight Matrix
gcdist <- ge.clst %>% 
  dplyr::mutate(longnum = sf::st_coordinates(geometry)[,1],
                latnum = sf::st_coordinates(geometry)[,2]) %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  sf::st_distance(x = .,
                  which = "Great Circle")

gcdist  <- apply(gcdist, 2, as.numeric) # drop units
gcdist.inv <- 1/gcdist
diag(gcdist.inv) <- 0

mp.clst$moranIclust <- purrr::map(mp.clst$data, 
                                  moran_mc_wrapper, 
                                  listw = gcdist.inv)

moranIclst <- mp.clst %>% 
  dplyr::select(c("plsmdmspec", "moranIclust")) 


moranIclst$Istatistic <- purrr::map(moranIclst$moranIclust, "statistic")
moranIclst$pvalue <- purrr::map(moranIclst$moranIclust, "p.value")
#......................
# save out
#......................
saveRDS(moranIclst, "results/MoranI_clust.RDS")
