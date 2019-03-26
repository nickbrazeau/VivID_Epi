#----------------------------------------------------------------------------------------------------
# SEASONALITY
#----------------------------------------------------------------------------------------------------

# adapted from OJWatson/magenta
# extracted from Cairns et al. 2012 PMC3621394 

calc_seasonality <- function(drc_admin_units_seasonal) {
  
  
  ## extract fourier covariates
  ssa0 <- drc_admin_units_seasonal$seasonal_a0
  ssa1 <- drc_admin_units_seasonal$seasonal_a1
  ssa2 <- drc_admin_units_seasonal$seasonal_a2
  ssa3 <- drc_admin_units_seasonal$seasonal_a3
  ssb1 <- drc_admin_units_seasonal$seasonal_b1
  ssb2 <- drc_admin_units_seasonal$seasonal_b2
  ssb3 <- drc_admin_units_seasonal$seasonal_b3
  theta_c <- rowMeans( drc_admin_units_seasonal[,grepl("seasonal_", colnames(drc_admin_units_seasonal))] )
    
  # define vector of times spanning one year
    tvec = 1:365
    
    # calculate Fourier series
    y <- sapply(tvec, function(x) {
      max((ssa0 + 
             ssa1 * cos(2 * pi * x / 365) + 
             ssa2 * cos(2 * 2 * pi * x / 365) + 
             ssa3 * cos(3 * 2 * pi * x / 365) + 
             ssb1 * sin(2 * pi * x / 365) + 
             ssb2 * sin(2 * 2 * pi * x / 365) + 
             ssb3 * sin(3 * 2 * pi * x / 365)) / theta_c, 0.001) # note the max 0.0001 ensures that scaling factor never goes below zero (this can happen in practice because we are only using the first few terms in an infinite series)
    })
  
    # now have 365 days of seasonality
  return(y)
}

# now add by month and admin level
# still adapted from OJWatson/magenta

####### 
#Issue with how regions are being linked
####### 


add_seasonal <- function(fourier_transform_365days, months, admin1){
    if(is.na(months)){
      return(NA)
    }
    months <- as.numeric(months)
    region <- names(attr(admin1,"labels"))[as.numeric(admin1)]
    ring <- c(12,1:11)
    mts <- ring[months]
    days <- lubridate::days_in_month(1:12) %>% cumsum
    days <- sapply(1:12,function(x) c(0,days[-12])[x]:days[x])
    return(mean(fourier_transform_365days[unlist(days[mts]),region]))
  }
  
  
  
# dt$seasonal_adm_month <- apply(dt,
#                                1,
#                                seasonal_add,   
#                                   fourier_transform_365days = calc_seasonality(drc_admin_units_seasonal),
#                                   months = dt$hv006,
#                                   admin1 = dt$hv024)
# 


#----------------------------------------------------------------------------------------------------
# Read Raster Temp and Precip
#----------------------------------------------------------------------------------------------------
readRasterBB <- function(rstfile, bb = bb){
  ret <- raster::raster(rstfile)
  ret <- raster::crop(x = ret, y = sf::as_Spatial(bb))
  ret <- raster::projectRaster(from = ret, to = ret,
                               crs = sf::st_crs("+proj=utm +zone=34 +datum=WGS84 +units=m")) # want units to be m
  
  return(ret)
  
}

extractDHSpts_from_raster <- function(hv001, raster, 
                                      geometry, urban_rura){
  # note raster is expected to have units of m
  
  geometry <- sf::as_Spatial(geometry)
  buffer = ifelse(urban_rura == "R", 10, 2)
  ret <- raster::extract(x = raster[[1]],
                         y =  geometry[[1]], 
                         buffer = buffer,
                         fun = mean,
                         sp = F) 
  return(ret)
}






                                 