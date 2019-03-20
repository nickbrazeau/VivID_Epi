
# CHECK IF TOOK ANY ANTIMALARIAL 2 WEEKS BEFORE COMES OUT

# read in PR



kds <- pr %>% 
  filter(hml16a < 60) %>% 
  left_join(x=., y=ge) %>% 
  dplyr::filter(latnum != 0 & longnum != 0) %>% 
  dplyr::filter(!is.na(latnum) & !is.na(longnum)) %>% 
  dplyr::mutate(hv005_wi = hv005/1e6,
                hml32_fctb = haven::as_factor(hml32),
                hml32_fctb = if_else(hml32_fctb %in% c("positive", "negative"), 
                                     hml32_fctb, factor(NA)),
                hml32_fctb =  forcats::fct_drop(hml32_fctb),
                hml32_numb = ifelse(hml32_fctb == "positive", 1, ifelse(hml32_fctb == "negative", 0, NA)),
                
                hml35_fctb = haven::as_factor(hml35),
                hml35_fctb = if_else(hml35_fctb %in% c("positive", "negative"), 
                                     hml35_fctb, factor(NA)),
                hml35_fctb =  forcats::fct_drop(hml35_fctb),
                hml35_numb = ifelse(hml35_fctb == "positive", 1, ifelse(hml35_fctb == "negative", 0, NA))
  )

options(survey.lonely.psu="certainty")
kdsrvy <- kds %>% srvyr::as_survey_design(ids = hv001, 
                                          strata = hv023, weights = hv005_wi)

rdtmicro <- kdsrvy %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(hv001) %>% 
  dplyr::summarise(
                   RDTprev = srvyr::survey_mean(hml35_numb, na.rm = T, vartype = c("se", "ci"), level = 0.95),
                   microprev = srvyr::survey_mean(hml32_numb, na.rm = T, vartype = c("se", "ci"), level = 0.95)
  ) %>% 
  mutate_if(is.numeric, round, 2)







