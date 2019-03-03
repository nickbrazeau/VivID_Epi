load("~/Documents/GitHub/VivID_Epi/data/vividmaps_small.rda")

dt$hv270_fctm <- forcats::fct_relevel(dt$hv270_fctm, c("poorest", "poorer", "middle", "richer", "richest"))
m <- xtabs(~dt$hv270_fctm + dt$pv18s + dt$hv025_fctb)


fisher.test(m[,,1])
fisher.test(m[,,2]) # has a one cell count


prop.table(m[,,1], margin = 1)
prop.table(m[,,2], margin = 1)

vars = c('hv270_fctm', 'hv025_fctb')

tableone::CreateTableOne(data=dt, 
                         strata = "pv18s", 
                         vars = vars,
                         includeNA = T)

tableone::CreateTableOne(data=dt, 
                         strata = c("hv025_fctb"),
                         vars = "hv270_fctm",
                         includeNA = T)
mantelhaen.test(m)



round(
  prop.table(
    ftable(dt[,c("hv025_fctb", "hv270_fctm", "pfldh")]),
    margin = 1),
3)

round(
  prop.table(
    ftable(dt[,c("hv025_fctb", "pfldh")]),
    margin = 1),
  3)


richpeopleinruralclusters <- dt[which(dt$hv025_fctb == "rural" & dt$hv270_fctm == "richest"), ]

richpeopleinruralclusterspospv <- dt[which(dt$hv025_fctb == "rural" & dt$hv270_fctm == "richest" & dt$pv18s == 1), ]

ggplot() + 
  geom_sf(data = DRCprov) +
  geom_point(data = richpeopleinruralclusters, aes(x=longnum, y=latnum)) +
  geom_point(data = richpeopleinruralclusterspospv, aes(x=longnum, y=latnum), color = "red") +
  vivid_theme +
  theme(axis.text = element_blank(),
        axis.line = element_blank())





chisq.test(dt$hv270_fctm, dt$pv18s)
chisq.test(dt$hv025_fctb, dt$pv18s)

chisq.test(dt$hv025_fctb, dt$hv270_fctm)


dt$hv270_fctm <- forcats::fct_relevel(dt$hv270_fctm, "poorest")


# model comparison

# urbanicity binary 
m1 <-  glm(pv18s ~ hv025_fctb,
              data = dt,
              family = binomial(link = 'logit'))
broom::tidy(m1, exponentiate=TRUE, conf.int=TRUE)


# wealth binary 
m2 <-  glm(pv18s ~ hv270_fctb,
           data = dt,
           family = binomial(link = 'logit'))
broom::tidy(m2, exponentiate=TRUE, conf.int=TRUE)

# wealth multi 
m3 <-  glm(pv18s ~ hv270_fctm,
           data = dt,
           family = binomial(link = 'logit'))
broom::tidy(m3, exponentiate=TRUE, conf.int=TRUE)

# urban confounding wealth b 
m4 <-  glm(pv18s ~ hv270_fctb + hv025_fctb,
           data = dt,
           family = binomial(link = 'logit'))
broom::tidy(m4, exponentiate=TRUE, conf.int=TRUE)


# wealth b confounding urban (same model)
m5 <-  glm(pv18s ~ hv025_fctb + hv270_fctb,
           data = dt,
           family = binomial(link = 'logit'))
broom::tidy(m5, exponentiate=TRUE, conf.int=TRUE)




# poor stratum, this is just the model above
c1 <- contrast::contrast(model,
         a=list(hv270_fctb = "poor_b", hv025_fctb = "rural"),
         b=list(hv270_fctb = "rich_b",  hv025_fctb = "rural"))

exp(c1$Contrast)
exp(c1$Lower)
exp(c1$Upper)

# rich stratum 
c2 <- contrast::contrast(model,
         a=list(hv270_fctb = "rich_b", hv025_fctb = "urban"),
         b=list(hv270_fctb = "poor_b",  hv025_fctb = "urban"))

exp(c2$Contrast)
exp(c2$Lower)
exp(c2$Upper)





### TYPE OF CLUSTER
round(
prop.table(
  ftable(dt[,c("hv025_fctb", "hv026_fctb", "pv18s_fctb")]),
  margin = 1),
3)
tableone::CreateTableOne(data=dt, 
                         strata = "pv18s_fctb", "hv026_fctm", includeNA = T)

dt %>% 
  group_by(hv026_fctm) %>% 
  summarise(clst = sum(!duplicated(hv001)))
  
broom::tidy(glm(pv18s ~ hv026_fctm + hv270_fctm,
           data = dt,
           family = binomial(link = 'logit')), 
           exponentiate=TRUE, conf.int=TRUE)

dt %>% 
    group_by(hv026_fctb) %>% 
    summarise(clst = sum(!duplicated(hv001)))

lrgcty <- dt[dt$hv026_fctm == "capital, large city", ]

ggplot() + 
  geom_sf(data = DRCprov) +
  geom_point(data = dt, aes(x=longnum, y=latnum, group = hv026_fctm, colour = hv026_fctm)) +
  vivid_theme +
  theme(axis.text = element_blank(),
        axis.line = element_blank())


