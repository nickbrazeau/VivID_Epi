#------------------------------------------
# Occupation
#------------------------------------------
# Going to dichotomize as suspected outdoor versus indoor
# Will code "others" as NA because cannot discern, only 2 obs

dt$hab717 <- ifelse(haven::as_factor(dt$hv104) == "female", dt$v717, dt$mv717)
table(dt$v717)
table(dt$mv717)
table(dt$hab717)

occupation <- readr::read_csv("~/Documents/GitHub/VivID_Epi/data/internal_datamap_files/pr_occupation_liftover.csv")
occupation$hab717 <- factor(occupation$hv717) # my mac being a pain

dt <- dt %>% 
  dplyr::mutate(hab717 = haven::as_factor(hab717), 
                hab717 =  forcats::fct_drop(hab717))
dt <- dt %>%
  left_join(x=., y=occupation, by="hab717") %>% 
  dplyr::mutate(hab717_fctb = factor(jobconditions),
                hab717_fctb = forcats::fct_relevel(hab717_fctb, "indoors"))

xtabs(~ dt$hab717 + dt$hab717_fctb, addNA = T)


