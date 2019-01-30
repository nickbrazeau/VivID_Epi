library(raster)
worldclim_prec = getData(name = "worldclim", var = "prec", res = 10)



ggplot() + geom_sf(data = DRCprov)

ggmap(drc_back_terrain) + 
  geom_sf(data = DRCprov, inherit.aes = F, alpha=0.99)



  geom_point(data = pv18sclust, aes(x = longnum, y = latnum, 
                                               fill = plsmd, colour = plsmd, size = plsmdn),
                        alpha = 0.8) +
    scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000") + 
    # scale_color_gradient2("Prevalence", low = "#0000FF", mid = "#FFEC00", high = "#FF0000", midpoint = quantile(data$plsmd[data$plsmd != 0], 0.75)) + 
    scale_size(guide = 'none') +  scale_fill_continuous(guide = 'none') +
    coord_sf(datum=NA)  # to get rid of gridlines
