#---------------------------------------------------------------------------------
# set up DRC borders
#---------------------------------------------------------------------------------
# https://github.com/ropensci/osmdata

bb <- osmdata::getbb("Democratic Republic of the Congo", 
                     featuretype = "country",
                     format_out = "sf_polygon")

