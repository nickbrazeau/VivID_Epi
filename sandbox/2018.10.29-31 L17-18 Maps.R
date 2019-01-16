#............................................
## Maps 1 - libraries & ecosystems ####
#............................................
# load("D:/User/Dropbox (Personal)/Education/Classes/18Fall_EPID799C_RforEpi/map_parts.rdata")
library(tidyverse)
library(broom)
library(sp)
library(rgdal) #for reading and writing shapefiles
library(rgeos) #for spatial tools like gCentroid
library(RColorBrewer) #For better, much better color ramps
library(tmap)
library(sf) #if you use the sf way, pretty much just need sf
library(ggthemes) # for theme_map
library(mapproj)

head(county_results) #formatter[formatter$variable=="cores",]

#.........................
# Get the county shapefile
#.........................
getwd()
if(!file.exists("./maps/nc counties.shp")){ #Run Once to get NC Counties.
  dir.create("./maps/")
  # ftp://ftp2.census.gov/geo/tiger/TIGER2016/ or
  # https://www.census.gov/cgi-bin/geo/shapefiles/index.php
  download.file("ftp://ftp2.census.gov/geo/tiger/TIGER2015/COUNTY/tl_2015_us_county.zip", "./maps/uscounties.zip")
  unzip("./maps/uscounties.zip", exdir="./maps")
  dir("./maps/")
  us_counties = readOGR(dsn="./maps", layer="tl_2015_us_county", stringsAsFactors = F)
  plot(us_counties)
  nc_counties = us_counties[us_counties$STATEFP == 37,] #look how slick this is.
  writeOGR(nc_counties, dsn="./maps", layer="nc counties", driver="ESRI Shapefile", overwrite_layer = T)
}

#.........................
# sp way
#.........................
nc_counties = readOGR(dsn="./maps", layer="nc counties", stringsAsFactors = F) # the sp way
class(nc_counties)
plot(nc_counties)
str(nc_counties, max.level = 2)
# ... or us tigris package. Also see acs package

#remember to merge with the spatial object, not just the data slot, to keep it from reordering!
nc_counties = merge(nc_counties, county_results, by.x="GEOID", by.y="FIPS")
head(nc_counties@data)
# Quick sp plots
spplot(nc_counties, "pct_pnc5", main="% Early PNC") # on spplot: https://edzer.github.io/sp/
spplot(nc_counties, "pct_preterm", main="% Preterm")
spplot(nc_counties, c("pct_pnc5", "pct_preterm"))
writeOGR(nc_counties, dsn="./maps", layer="nc counties w data", driver="ESRI Shapefile", overwrite_layer = T)

sum(nc_counties@data$n)
sum(nc_counties[nc_counties$NAME %in% c("Orange", "Durham", "Chatham"),]$n)

### A little prettier
nc_counties$pct_earlyPNC_f = cut(nc_counties$pct_pnc5, 7)
spplot(nc_counties, "pct_earlyPNC_f",
       col.regions=brewer.pal(7,"RdYlGn"),
       main="% of Births with Early Start of Prenatal Care (<5mo)")

### Using base plot()
display.brewer.all()
nc_counties$pct_preterm_colors = as.character(cut(nc_counties$pct_preterm, 10, labels = brewer.pal(10,"RdYlGn")))
plot(nc_counties, co=nc_counties$pct_preterm_colors, main="% Preterm birth")
county_centroids = gCentroid(nc_counties, byid = T) # aaand here's why R's great as a GIS.
text(county_centroids@coords, labels=as.character(nc_counties$NAME), cex=0.4)

### Using ggplot() OLD style - nice, but slow
names(nc_counties@data)
nc_counties_fort = fortify(nc_counties, region = "GEOID") # This is very old fashioned
str(nc_counties_fort)
nc_counties_fort = merge(nc_counties_fort, nc_counties, by.x="id", by.y="GEOID")
head(nc_counties_fort)
cut(nc_counties$pct_pnc5, 5)
pretty(nc_counties$pct_pnc5) #funny name
nc_counties_fort$pct_pnc5_f = cut(nc_counties_fort$pct_pnc5, breaks=pretty(nc_counties_fort$pct_pnc5))
table(nc_counties_fort$pct_pnc5_f)
#grep to make , a -_ move to top, like
#http://time_com/4394141/zika-abortion-usa/
mymap =
  ggplot(nc_counties_fort, aes(x=long, y=lat, group=group, fill=pct_pnc5_f)) +
  geom_polygon(color="white")+coord_map()+
  labs(title="Early prenatal care by county in North Carolina, 2012",
       subtitle="This ggplot is gorgeous, holy cow",
       x="", y="")+
  scale_fill_brewer(name="% Early Prenatal Care", type="seq", palette ="Greens")+
  #scale_fill_gradient(name="% Early Prenatal Care", low = "white", high = "dark green")+
  theme_map()+theme(legend.position="right")
#theme_map()+theme(legend_position="bottom")+guides(fill=guide_legend(nrow=1))#right")
mymap
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
mymap+theme_minimal() #also a nice one / clean.

# tmap
tm_shape(nc_counties)+ # skip the whole fortify thing.
  tm_fill("pct_pnc5", title="% Early PNC")+
  tm_borders()+ #could add more layers with tm_shape and then tm_whatever...
  tm_layout(main.title="tmap: % Early Prenatal Care by County")
# ^ Note ggplot like syntax.

# sf - the current best practice
library(sf)
nc_counties_sf = st_read("./maps/nc counties.shp", stringsAsFactors = F)
str(nc_counties_sf) # a dataframe
plot(nc_counties_sf)

nc_counties_sf %>% # dplyr the df right into ggplot!
  mutate(GEOID = as.numeric(GEOID)) %>%
  left_join(county_results, by=c("GEOID"="FIPS")) %>%
  ggplot(aes(fill = pct_pnc5))+
  geom_sf()+ # default is that geometry "aesthetic" is in geometry column. Which it is!
  theme_minimal() + # alternatives here
  labs(title="sf and ggplot - perhaps best for now: % Early PNC ")
# make_EPSG() # or lookup

nc_counties_sf %>% # dplyr the df right into ggplot!
  mutate(GEOID = as.numeric(GEOID)) %>%
  left_join(county_results, by=c("GEOID"="FIPS")) %>%
  mutate(preterm_quintile = cut(pct_preterm, quantile(pct_preterm, seq(0,1,.2)))) %>%
  group_by(preterm_quintile) %>% # how cool is this.
  summarise(pct_pnc5 = mean(pct_pnc5, na.rm=T)) %>%
  ggplot(aes(fill = pct_pnc5))+
  geom_sf()
#................................

#................................
# Maps 2 - spatial analysis ####
#................................
# Demos: transform, gTouches,
proj4string(nc_counties) #or nc_counties@proj4string
nc_stateplane_proj = "+proj=lcc +lat_1=36.16666666666666 +lat_2=34.33333333333334 +lat_0=33.75 +lon_0=-79 +x_0=609601.2192024384 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs"
nc_counties_stateplane = spTransform(nc_counties, nc_stateplane_proj)

# centroids
county_cents = gCentroid(nc_counties_stateplane, byid=T) #just the spatial object
county_cents_sp = SpatialPointsDataFrame(county_cents, data=nc_counties@data) #need to add back the data

# buffers
plot(nc_counties_stateplane); plot(county_cents, add=T)
plot(county_cents_sp[county_cents_sp$NAME == "Orange",], col="blue", add=T, pch="0")
plot(gBuffer(county_cents, width = 10*5280, byid = T), add=T) #10 mi buffers

# spatial relationships - touching?
plot(nc_counties[county_cents_sp$NAME == "Orange",])
orange_neighbors = gTouches(nc_counties[county_cents_sp$NAME == "Orange",], nc_counties, byid = T)
str(orange_neighbors) # note it's a matrix
plot(nc_counties); plot(nc_counties[orange_neighbors[,1],], col="blue", add=T)

# A question of hospitals
# http://data.nconemap.gov/geoportal/rest/find/document?searchText=Hospitals&f=searchpage&f=searchpage
med_facilities = read_sf("./maps/medfacs.shp", stringsAsFactors = F) # the sf way
plot(nc_counties); plot(med_facilities %>% st_geometry(), add=T)
nc_counties_sf = nc_counties %>% st_as_sf()
st_crs(nc_counties_sf); st_crs(med_facilities) # oops! projection mismatch!

med_facilities = med_facilities %>% st_transform(crs = nc_stateplane_proj)
nc_counties_sf = nc_counties_sf %>% st_transform(crs = nc_stateplane_proj)
plot(nc_counties_sf %>% st_geometry()); plot(med_facilities %>% st_geometry(), add=T, col="blue", pch=".")

#what centroids are 50m from a hospital?
head(med_facilities)
table(med_facilities$STYPE)
hospitals = med_facilities[med_facilities$STYPE == "Hospital",] # nice!

# don't have birth locations - let's use centroids as a proxy
dist_matrix = st_distance(county_centroids %>% spTransform(nc_stateplane_proj) %>% st_as_sf, hospitals)
dim(dist_matrix)
dist_matrix[1,]
within_5mi = function(x){any(x<=5*5280)}
apply(dist_matrix, 1, within_5mi) # can colbind this back together.

dist_list = st_is_within_distance(county_centroids %>% spTransform(nc_stateplane_proj) %>% st_as_sf,
                                  hospitals, 5*5280)
map_lgl(dist_list, ~ length(.x)>0)


# spatial over... so many ways
hold = st_within(hospitals, nc_counties_sf)
str(hold) # what do we know that can work with lists?
unlist(hold) # here's one way!
map_int(hold, ~.x[[1]]) # here's another flexible way.

save(list = c("nc_counties", "nc_counties_sf", "county_results", "med_facilities", "nc_stateplane_proj"), file = "map_parts.rdata")
#.................................................