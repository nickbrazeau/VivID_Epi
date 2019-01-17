#----------------------------------------------------------------------------------------------------
# Purpose of this script is to investigate if the pv and pf data are correlated
# considering both space and other epi variables
#----------------------------------------------------------------------------------------------------



# Comparison of _Pv_ and _Pf_ Risk Factors
## Age {.tabset .tabset-fade .tabset-pill}
```{r}

pvage <- dt %>% 
  dplyr::mutate(count = 1) %>% 
  srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
  dplyr::group_by(hv105) %>% 
  dplyr::summarise(plsmdn = srvyr::survey_total(count), 
                   plsmd = srvyr::survey_mean(pv18s, na.rm = T, vartype = c("se", "ci"), level = 0.95))
pvage$species <- "vivax"

pfage <- dt %>% 
  dplyr::mutate(count = 1) %>% 
  srvyr::as_survey_design(ids = hv001, weights = hiv05_cont) %>% 
  dplyr::group_by(hv105) %>% 
  dplyr::summarise(plsmdn = srvyr::survey_total(count), 
                   plsmd = srvyr::survey_mean(pfldh, na.rm = T, vartype = c("se", "ci"), level = 0.95))

pfage$species <- "falciparum"

plotObj <- dplyr::bind_rows(pvage, pfage) %>% 
  ggplot() + 
  geom_line(aes(x=hv105, y=plsmd, color = species)) + 
  geom_ribbon(aes(x=hv105, ymin=plsmd_low, ymax=plsmd_upp, fill = species), alpha=0.5) +
  geom_point(aes(x=hv105, y=plsmd, size=plsmdn, color = species), alpha=0.5, show.legend=F) +
  scale_fill_manual("Malara Species", values = c("#a50026", "#313695"), labels = c("Pfalciparum", "Pvivax")) +
  scale_colour_manual(values = c("#a50026", "#313695"), guide =F) +
  ggtitle("Malaria Prevalence by Age") +
  xlab("Age") + ylab("Malaria Prevalence") + 
  vivid_theme

```
## 

