## .................................................................................
## Purpose: Ape Permutation testing
##
## Notes: 
## .................................................................................
set.seed(48)
library(tidyverse)
library(sf)
source("R/00-functions_basic.R")

#............................................................
# Import Data
#............................................................
load("data/map_bases/vivid_maps_bases.rda")
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
DRCprov <- readRDS("data/map_bases/vivid_DRCprov.rds")
ge <- readRDS(file = "data/raw_data/dhsdata/VivIDge.RDS")

# need to add back in geometry to dt
gegeoms <- ge %>% 
  dplyr::select(c("hv001", "geometry"))
dt <- dplyr::left_join(dt, gegeoms, by = "hv001")

dtsrvy <- makecd2013survey(survey = dt)
ape <- readRDS("data/derived_data/drc_ape.rds")
# sanity check
sf::st_crs(ape)

#............................................................
# set up permuations
#...........................................................
# how many clusters do apes overlap?
ape.range <- sf::st_union(ape$geometry)
ape.clusters <- sf::st_intersection(ge, ape.range)

# how many clusters do chimp and gorillas overlap?
chimpgor.range <- ape %>% 
  dplyr::filter(species != "Pan paniscus") %>% 
  dplyr::select(geometry) %>% 
  sf::st_union(.)
chimpgor.clusters <- sf::st_intersection(ge, chimpgor.range)


#............................................................
# run permutations
#...........................................................
run_ape_permutations <- function(dtsrvy18s, n.apeclsts){
  
  clsts <- dtsrvy18s %>% 
    tibble::as_tibble(.) %>% 
    dplyr::select("hv001") %>% 
    dplyr::filter(!duplicated(.)) %>% 
    unlist(.)
  
  apeclsts <- sample(clsts, size = n.apeclsts, replace = F)
  
  ape.ovrlp.pos <- dtsrvy18s %>% 
    dplyr::filter(hv001 %in% apeclsts) %>% 
    dplyr::summarise(pv18s = srvyr::survey_mean(pv18s)) %>% 
    dplyr::select(-c(dplyr::ends_with("_se")))
  
  ret <- unname(as.vector(ape.ovrlp.pos))
  
  return(ret)
  
}

sims <- 1e4


#......................
# for all apes
#......................
ape.range.iters <- data.frame(iter = 1:sims)

ape.range.iters$ape_overlap_prev <- purrr::map(ape.range.iters$iter, function(x){
  return( run_ape_permutations(dtsrvy18s = dtsrvy, n.apeclsts = length(ape.clusters)) )
}) 

ape.range.iters$ape_overlap_prev <- unlist(ape.range.iters$ape_overlap_prev )

# save out
dir.create("analyses/08-nonparam_riskfacts/perm_rets/", recursive = T)
saveRDS(ape.range.iters, "analyses/08-nonparam_riskfacts/perm_rets/all_ape_range_iters.RDS")


#......................
# for all chimps and gorillas 
#......................
chimpgor.range.iters <- data.frame(iter = 1:sims)

chimpgor.range.iters$chimpgor_overlap_prev <- purrr::map(chimpgor.range.iters$iter, function(x){
  return( run_ape_permutations(dtsrvy18s = dtsrvy, n.apeclsts = length(chimpgor.clusters)) )
}) 

chimpgor.range.iters$chimpgor_overlap_prev <- unlist(chimpgor.range.iters$chimpgor_overlap_prev )

# save out
saveRDS(chimpgor.range.iters, "analyses/08-nonparam_riskfacts/perm_rets/chimp_gor_range_iters.RDS")


