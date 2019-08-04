#----------------------------------------------------------------------------------------------------
# Purpose of this script is to collect power calculations
#----------------------------------------------------------------------------------------------------

#........................... 
# Read in power results
#........................... 
poweriters.files <- list.files(path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/09-Power/_rslurm_powercalcs",
                              pattern = ".RDS", full.names = T) 
# reorder file paths correctly
poweriters.files <- tibble::tibble(files = poweriters.files) %>% 
  dplyr::mutate(filename = basename(as.character( files )),
                order = stringr::str_extract(filename, "[0-9]+"),
                order = as.numeric(order)) %>%
  dplyr::filter(!is.na(order)) %>% 
  dplyr::arrange(order) %>% 
  dplyr::select(files) %>% unlist(.)

poweriters.ret <- purrr::map(poweriters.files, function(x){
                                                            ret <- readRDS(x)
                                                            ret.df <- dplyr::bind_rows(ret)
                                                            return(ret.df)
                                                            }) %>% 
  dplyr::bind_rows()

#........................... 
# Read in power params
#........................... 
poweriters.params <- readRDS("~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/09-Power/_rslurm_powercalcs/params.RDS")
poweriters.params.ret <- dplyr::bind_cols(poweriters.params, poweriters.ret)


#........................... 
# Get Summary
#........................... 

poweriters.powercalc <- poweriters.params.ret %>% 
  dplyr::group_by(n, p, exp_prob, p0, RR) %>% 
  dplyr::summarise(
    power = sum(p1 < 0.05)
  )



jpeg(filename = "results/figures/RR_glm_posthoc_powercalc.jpg", width = 8, height = 6, res = 500, units = "in")

poweriters.powercalc %>%
  dplyr::mutate(beta = 1-power,
                expprob_f = factor(exp_prob)) %>%
  ggplot() +
  geom_line(aes(x=beta, y=RR, color = expprob_f), method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
  geom_point(aes(x=beta, y=RR, color = expprob_f), stat="identity") +
  geom_vline(aes(xintercept=0.2), colour="#de2d26", linetype="dashed") +
  ggtitle(label="Simulated Risk Ratio versus \n Type II Error (Complement of Power)") +
  xlab("Type II Error") + ylab("Odds Ratio") +
  scale_color_manual("Exposure Probability", values = c("#8214A0", "#005AC8", "#006E82")) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.text = element_text(hjust = 0.5, size=13))



graphics.off()


