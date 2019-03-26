library(GSODR)
library(future)
library(PrevMap)
#----------------------------------------------------------------------------------------------------
# Extract station information from GSODR
#----------------------------------------------------------------------------------------------------
# https://ropensci.org/blog/2017/04/04/gsodr/
findneareststations <- function(input = dt, urbdist = 50, rurdist = 250, latestdate = 20141231){
  cd2013latlong <- data.frame(
    hv001 = dt$hv001,
    LAT = dt$latnum,
    LON = dt$longnum,
    ruralurban = haven::as_factor(dt$hv025)
  )
  cd2013latlong <- cd2013latlong[!duplicated(cd2013latlong),]
  cd2013latlong <- cd2013latlong %>% 
    dplyr::mutate(distance = ifelse(!is.na(ruralurban) & ruralurban == "urban",
                                    urbdist, 
                                    ifelse(!is.na(ruralurban) & ruralurban == "rural",
                                           rurdist, NA
                                    )))
  
  # Fetch station list from NCEI (thanks ROpenSci)
  # Need this metadata to see when last station data was
  station_meta <- read_csv(
    "ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv",
    col_types = "ccccccddddd",
    col_names = c("USAF", "WBAN", "STN_NAME", "CTRY", "STATE", "CALL", "LAT",
                  "LON", "ELEV_M", "BEGIN", "END"), skip = 1)
  station_meta$STNID <- as.character(paste(station_meta$USAF,
                                           station_meta$WBAN,
                                           sep = "-")) 
  
  nrstat <- cd2013latlong %>% 
    dplyr::select(-c("hv001", "ruralurban")) %>% 
    purrr::pmap(., GSODR::nearest_stations) # using haversine distance -- good
  
  nrstat <- lapply(nrstat, function(x){
    stout <- x[ x %in% station_meta$STNID[station_meta$END >= latestdate] ]
    stout <- x[x != "641260-99999"] # hard code this station that is missing data despite end date being ok
    return(stout)
    
  })
  
  
  ret <- data.frame(
    hv001 = cd2013latlong$hv001,
    ruralurban = cd2013latlong$ruralurban
  ) 
  ret$nrstat <- nrstat
  if(any( sapply(ret$nrstat, function(x){return(length(x) == 0)}) )){
    stop(paste("Cluster", paste0(ret$hv001[which(sapply(ret$nrstat, function(x){return(length(x) == 0)}))], ","),
               "a", ret$ruralurban[which(sapply(ret$nrstat, function(x){return(length(x) == 0)}))], 
               "cluster, did not find a station. Increase search distance"))
  }
  
  out <- lapply(1:nrow(ret), function(x){
    df <- data.frame(
      hv001 = rep(ret$hv001[x], times = length(ret$nrstat[[x]])),
      STNID = ret$nrstat[[x]]
    )
    return(df)
  }) %>% 
    dplyr::bind_rows(.) 
  
  # end of function
  return(out)
}

pull_station <- function(stations, years){
  future::plan("multisession")
  gsod <- get_GSOD(years = years,
                   station = stations)
  return(gsod)
}


nrststat <- findneareststations(input = dt, urbdist = 50, rurdist = 250, latestdate = 20141231)
gsod <- pull_station(stations = nrststat$STNID[!duplicated(nrststat$STNID)], 
                     years = 2013:2014)
tempprecip <- dplyr::left_join(x = nrststat, y = gsod)

tempprecip <- tempprecip %>% 
  dplyr::mutate(YEARMODA = lubridate::ymd(YEARMODA)) %>% 
  dplyr::group_by(hv001, YEARMODA) %>% 
  dplyr::summarise(meanDailyTemp = mean(TEMP),
                   meanDailyPrecip = mean(PRCP))





