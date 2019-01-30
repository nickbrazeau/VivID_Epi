library(shp2graph) #Shp2graph: Tools to Convert a Spatial Network into an Igraph Graph in R -- manuscript
library(tidygraph)
library(ggraph)
load("~/Documents/GitHub/VivID_Epi/data/osm_roads.rda")

bb <- getbb("Kinshasa")
bb <- st_bbox(st_sf(bb = 1:2, geom = st_sfc( st_point(bb[,1]), st_point(bb[,2]) ) ))

init <- clstrs$data[[1]]
init <- sf::st_as_sf(prev)

rds <- sf::st_crop(rds, bb)
init <- sf::st_crop(init, bb)

pts <- as.matrix(init[,c("longnum", "latnum")])[,1:2] # odd behavior bc it keep geometry ... not like a dataframe

t <- shp2graph::points2network(ntdata = sf::as_Spatial(road),
                          pointsxy = points,
                          approach = 1) #https://rdrr.io/cran/shp2graph/man/points2network.html
t.cnt <- nt.connect(sf::as_Spatial(road))
plot(t.cnt)
u <- shp2graph::nel2igraph(nodelist = t[[1]], edgelist = t[[2]])
plot(u)

igraph::shortest_paths(u, from = t)
shortest.paths(u, v=V(u), to=V(u))

v <- tidygraph::as_tbl_graph(u)
