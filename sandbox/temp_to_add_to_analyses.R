load("~/Documents/GitHub/VivID_Epi/data/vividepi_recode.rda")
load("~/Documents/GitHub/VivID_Epi/data/04-basic_mapping_data.rda")
load("~/Documents/GitHub/VivID_Epi/data/osm_roads.rda")
load("~/Documents/GitHub/VivID_Epi/data/vividmaps_small.rda")

#..........................................
# setup
#..........................................
clustgeom <- dt[!duplicated(dt$hv001), c("hv001", "latnum", "longnum", "geometry")]
clstrs <- mp %>% 
  filter(maplvl == "hv001")
clstrs$data <- map(clstrs$data, function(x){
  return( dplyr::inner_join(x, clustgeom, by = "hv001") )
})

#..........................................
# mk dataframe
#..........................................
pvclst <- clstrs$data[[2]]
pfclst <- clstrs$data[[1]]
# looking at wealth and urbanicity
# going to treat hv270 wealth as ordinal even though it is a categorical variable (weakly ordinal argument)
wlthclst <- dt %>% 
  dplyr::group_by(hv001, hv026_fctm) %>% 
  dplyr::summarise(wlth = mean(hv270), longnum = mean(longnum), latnum=mean(latnum))

# joins
pvclst <- inner_join(pvclst, wlthclst, by = c("longnum", "latnum"))
pfclst <- inner_join(pfclst, wlthclst, by = c("longnum", "latnum"))






jpeg(filename = "~/Desktop/pv_wealth_urban.jpg", height = 8, width = 11, units = "in", res = 400)
ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = pvclst, aes(x=longnum, y=latnum, fill=wlth, shape = hv025_fctb, size = plsmdprev), 
             color = "grey20", alpha = 0.6) +
  viridis::scale_fill_viridis("Wealth Dist", option = "plasma") +
  scale_shape_manual("Rural v. Urban", values = c(21, 24)) +
  ggtitle("Pv Prev, Wealth, & Urbanicity by Cluster") +
  vivid_theme +
  theme(axis.text = element_blank(),
        axis.line = element_blank(), 
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 11, angle = 0),
        axis.title = element_blank())
graphics.off()




jpeg(filename = "~/Desktop/pf_wealth_urban.jpg", height = 8, width = 11, units = "in", res = 400)
ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = pfclst, aes(x=longnum, y=latnum, fill=wlth, shape = hv025_fctb, size = plsmdprev), 
             color = "grey20", alpha = 0.6) +
  viridis::scale_fill_viridis("Wealth Dist", option = "plasma") +
  scale_shape_manual("Rural v. Urban", values = c(21, 24)) +
  ggtitle("Pf Prev, Wealth, & Urbanicity by Cluster") +
  vivid_theme +
  theme(axis.text = element_blank(),
        axis.line = element_blank(), 
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 11, angle = 0),
        axis.title = element_blank())
graphics.off()




jpeg(filename = "~/Desktop/pv_wealth_city.jpg", height = 8, width = 11, units = "in", res = 400)
ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = pvclst, aes(x=longnum, y=latnum, fill=wlth, shape = hv026_fctm, size = plsmdprev), 
             color = "grey20", alpha = 0.6) +
  viridis::scale_fill_viridis("Wealth Dist", option = "plasma") +
  ggtitle("Pf Prev, Wealth, & Municipal-type by Cluster") +
  vivid_theme + 
  theme(axis.text = element_blank(),
        axis.line = element_blank(), 
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 11, angle = 0),
        axis.title = element_blank())
graphics.off()

