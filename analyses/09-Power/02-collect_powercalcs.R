#----------------------------------------------------------------------------------------------------
# Purpose of this script is to collect power calculations
#----------------------------------------------------------------------------------------------------

#........................... 
# Read in power results
#........................... 
load("results/powercalcs.rda")
poweriters.ret <- poweriters.ret %>% 
  dplyr::bind_rows()
poweriters.params.ret <- dplyr::bind_cols(poweriters.paramsdf, poweriters.ret)
#........................... 
# Get Summary
#........................... 

poweriters.powercalc <- poweriters.params.ret %>% 
  dplyr::group_by(n, p, exp_prob, p0, OR) %>% 
  dplyr::summarise(
    power = mean(p1 < 0.05)
  )

library(plotly)
df <- poweriters.powercalc %>%
  dplyr::mutate(beta = 1-power,
                expprob_f = factor(exp_prob))

jpeg(filename = "results/figures/OR_glm_posthoc_powercalc.jpg", width = 11, height = 8, res = 250, units = "in")


ggplot(df, aes(x=beta, y=OR, color = expprob_f)) +
#  stat_smooth(method = 'nls', formula = y ~ a * exp(x) + b, se = FALSE, start = list(a=-1, b=-1)) +
  geom_point() +
  geom_vline(aes(xintercept=0.2), colour="#de2d26", linetype="dashed") +
  ggtitle(label="Simulated Risk Ratio versus \n Type II Error (Complement of Power)") +
  xlab("Type II Error") + ylab("Odds Ratio") +
  scale_color_manual("Exposure Probability", values = c("#8214A0", "#005AC8", "#006E82")) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.text = element_text(hjust = 0.5, size=13))



graphics.off()


