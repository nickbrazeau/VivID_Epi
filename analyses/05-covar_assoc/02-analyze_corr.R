



#......................
# Parametric Approach -- will account for weights here since outcome dependent
#......................
dtsrvy <- makecd2013survey(survey = dt)
# VIF
eq <- as.formula(paste0("pv18s~", paste(rskfctr, collapse = "+")))
model.sat <- survey::svyglm(eq,
                            design = dtsrvy,
                            family = quasibinomial(link="logit"))
summary(model.sat)
vifs <- car::vif(model.sat)
summary(vifs)



#------------------------------------------------------------------------------------------
# Analyze for Missing Data
#-------------------------------------------------------------------------------------------
dtsrvy %>% 
  dplyr::select(rskfctr) %>% 
  as.data.frame(.) %>% 
  mice::md.pattern(., plot=T, rotate.names = T)

# Are my missing cases fundamentally different?
dt.cmpl <- dt %>% 
  dplyr::select(c("hv001", "hv023", "hv005_wi", rskfctr)) %>% 
  dplyr::mutate(
    ismissing = ifelse(rowSums(is.na(.)) == 0, "no", "yes"),
    ismissing = factor(ismissing))

sf::st_geometry(dt.cmpl) <- NULL
dt.cmpl.srvy <- makecd2013survey(dt.cmpl)


misstbl1 <- tableone::svyCreateTableOne(
  data = dt.cmpl.srvy,
  strata = "ismissing",
  vars = rskfctr,
  includeNA = T,
  smd = T,
  test = F)

print(misstbl1, smd = TRUE)

# are missing samples spread out?
dt.cmpl %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::filter(ismissing == "yes") %>% 
  dplyr::summarise(
    n = n()
  ) %>% 
  DT::datatable()


#................................
# Note to any reading this code: 
# This data does not appear to be MCAR... there do appear to be some biases (more or less, poor, rural people with lots of household members)
# However, given that this 0.44% of the unweighted data and 0.43% of the weighted data, we are well within the range of the rule of thumb
# that missing data less than 5-10% (even if it is MNAR) is largely inconsequential (see Bennett 2001). As a result, I feel comfortable subsetting to complete
# cases from here on out without worrying about the effect of any MNAR bias that may or may not be present and would affect my results
#................................

#......................
# Results & Out
#......................
saveRDS(corr.grid, file = "data/derived_data/covariate_correlation_data.rds")








