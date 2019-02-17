
# looking at wealth and urbanicity
# going to treat hv270 wealth as ordinal even though it is a categorical variable (weakly ordinal argument)
wlthclst <- dt %>% 
  dplyr::group_by(hv001, hv025_fctb) %>% 
  dplyr::summarise(wlth = mean(hv270), longnum = mean(longnum), latnum=mean(latnum))


jpeg(filename = "~/Desktop/wealth_urban.jpg", height = 8, width = 11, units = "in", res = 400)
ggplot() +
  geom_sf(data = DRCprov) +
  geom_point(data = wlthclst, aes(x=longnum, y=latnum, fill=wlth, shape = hv025_fctb), 
             color = "black", alpha = 0.6, size = 2) +
  viridis::scale_fill_viridis("Wealth Distribution", option = "plasma") +
  scale_shape_manual("Rural v. Urban", values = c(21, 24)) +
  ggtitle("Wealth versus Urbanicity by Cluster") +
  vivid_theme +
  theme(axis.text = element_blank(),
        axis.line = element_blank(), 
        legend.position = "right",
        legend.text = element_text(face = "bold", size = 11, angle = 0),
        axis.title = element_blank())
graphics.off()