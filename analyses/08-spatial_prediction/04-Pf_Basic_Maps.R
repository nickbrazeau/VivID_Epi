#----------------------------------------------------------------------------------------------------
# Purpose of this script is to create basic maps as sanity checks
# Will refer to Molly's paper for further Pf-mapping 
#----------------------------------------------------------------------------------------------------

# depending on what SRM says don't do this in a separeate script -- this is what pmap is for

#......................
# Import Data
#......................
mp <- readRDS("data/derived_data/basic_cluster_mapping_data.rds")


#......................
# Subset to Pfal
#......................
pf.prov.weighted <- mp %>% 
  dplyr::filter(plsmdmspec == "pfldh" & maplvl == "adm1name") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pf.prov.weighted <- sf::st_as_sf(pf.prov.weighted)
# need to keep integers, so will round
pf.prov.weighted <- pf.prov.weighted %>% 
  dplyr::mutate_if(is.numeric, round, 0)

pf.prov.weighted.nosf <- pf.prov.weighted
sf::st_geometry(pf.prov.weighted.nosf) <- NULL

pf.clust.weighted <- mp %>% 
  dplyr::filter(plsmdmspec == "pfldh" & maplvl == "hv001") %>% 
  dplyr::select(data) %>% 
  tidyr::unnest()
# vectors have destroyed spatial class, need to remake
pf.clust.weighted <- sf::st_as_sf(pf.clust.weighted)
# need to keep integers, so will round
pf.clust.weighted <- pf.clust.weighted %>% 
  dplyr::mutate_if(is.numeric, round, 0)
# drop sf obj for vectors
pf.clust.weighted.nosf <- pf.clust.weighted
sf::st_geometry(pf.clust.weighted.nosf) <- NULL


#-------------------------------------------------------------------------
# Conditional Autoregressive Spatial Model 
#-------------------------------------------------------------------------
#......................
# Make Adjacency Matrix for Pvf
#......................
# https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
W.nb <- spdep::poly2nb(sf::as_Spatial(pf.prov.weighted), row.names = pf.prov.weighted$adm1name)
W <- spdep::nb2mat(W.nb, style = "B") # binary weights taking values zero or one (only one is recorded)

pf.prov.ret <- CARBayes::S.CARleroux(formula = "plsmdn ~ 1", 
                                     W = W,
                                     data = pf.prov.weighted.nosf,
                                     trials = pf.prov.weighted.nosf$n, 
                                     prior.var.beta = 5e4,
                                     burnin = 1e3,
                                     n.sample = 1e5,
                                     family = "binomial")

#........................
# Moran's I for Pf Prov
#........................
pf.prov.moranI <- spdep::moran.mc(mp$data[[1]]$plsmdprev,
                                  listw = spdep::mat2listw(W),
                                  alternative = "greater",
                                  nsim = 1e5)




#----------------------------------------------------------------------------------------------------
# PrevMap & The Gaussian Process Model 
#----------------------------------------------------------------------------------------------------
# Following Giorgi & Diggle 2017 -- https://www.jstatsoft.org/article/view/v078i08
#......................
# EDA - Variogram
#......................
eda <- pf.clust.weighted.nosf %>% 
  dplyr::select(c("hv001", "plsmdn", "n", "plsmdprev", "longnum", "latnum"))

# transform count of "successes" to logit space
eda$plsmdnlogit <- log( (eda$plsmdn + 0.5)/(eda$n - eda$plsmdn + 0.5) ) # 0.5 as tolerance for 0s

# look at kappa smooth parameter for matern covariance
profile.kappa <- PrevMap::shape.matern(plsmdnlogit ~ 1,
                                       coords = ~ longnum + latnum,
                                       data = eda,
                                       set.kappa = seq(1e-3, 3, length = 16),
                                       start.par = c(0.2, 0.05), # starting values for the scale parameter phi and the relative variance of the nugget effect nu2
                                       coverage = NULL # CI coverage
)
 

# fit a voariogram 
coords <- eda[,c("longnum", "latnum")]
vari <- geoR::variog(coords = coords, data = eda$plsmdnlogit)
vari.fit <- geoR::variofit(vari, 
                           ini.cov.pars = c(3, 1),
                           cov.model = "matern",
                           fix.kappa = F)

plot(vari$v ~ vari$u, xlab = "distance", ylab = "semivariance")
plot(vari, main = substitute(paste("Semivariogram for ", italic("P. falciparum"))))
# really not an impressive amount of spatial autocorrelation




#........................
# Moran's I for Cluster Level
#........................
gc <- mp %>% 
  dplyr::filter(maplvl == "hv001" & plsmdmspec == "pfldh") %>% 
  tidyr::unnest() %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distHaversine) 

pfclust <- mp$data[[4]]
sf::st_geometry(pfclust) <- NULL
pfclust.dist <- pfclust %>% 
  dplyr::select(c("longnum", "latnum")) %>% 
  geosphere::distm(x =., fun = geosphere::distHaversine) 

pfclust.dist.inv <- 1/pfclust.dist
diag(pfclust.dist.inv) <- 0

mp.moranI.ret.gaus <- spdep::moran.mc(pfclust$plsmdprev,
                                      # temp to use MP,
                                      listw = spdep::mat2listw(pfclust.dist.inv),
                                      alternative = "greater",
                                      nsim = 1e5)












