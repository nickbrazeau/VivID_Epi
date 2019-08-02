


powersum <- data.frame(RR=numeric(), power=numeric())

for(k in 1:length(p0sim)){
  
  p0 <- p0sim[k]
  
  # running individual p0
  i <- 0
  ret <- data.frame(RR=numeric(), p=numeric())
  while(i<1000){
    temp <- powercalculator.glmRR(n=15879, exp_prob=0.5, p=0.03, p0=p0)
    ret <- rbind(ret,temp)
    i <- i +1
  }
  
  # saving results
  tempsum <- data.frame(RR=mean(ret$RR), power=sum(ret$p<0.05)/1000)
  powersum <- rbind(powersum, tempsum)
  
}







jpeg(filename = "results/figures/RR_glm_posthoc_exposure50perc_powercalc.jpg", width = 8, height = 6, res = 500, units = "in")

powersum %>%
  dplyr::mutate(beta = 1-power) %>%
  ggplot() +
  geom_point(aes(x=beta, y=RR), stat="identity") +
  geom_vline(aes(xintercept=0.2), colour="#de2d26", linetype="dashed") +
  ggtitle(label="Simulated Risk Ratio versus \n Type II Error (Complement of Power)") +
  xlab("Type II Error") + ylab("Odds Ratio") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.text = element_text(hjust = 0.5, size=13))



graphics.off()


# #----------
# # LOD from cedar
# lod <- read.csv(file="~/Desktop/pvlod.csv", header = T)
# levels(factor(lod$Sample.Name))
# decoder <- data.frame(Sample.Name=rev(c("Pv_0.00000003325", "Pv_0.0000000625", "Pv_0.000000125", "Pv_0.00000025",
#                          "Pv_0.0000005", "Pv_0.000001" )),
#            parasites=c("45.50",
#                        "22.75",
#                        "11.38",
#                        "5.69",
#                        "2.84",
#                        "1.42"))
# 
# 
# lod <- left_join(lod, decoder, by="Sample.Name")
# jpeg(filename = "~/Downloads/lod.jpg", width = 8, height = 6, res = 700, units = "in")
# lod %>%
#   dplyr::filter(CT != "Undetermined") %>%
#   dplyr::mutate(CTnum = as.numeric(as.character(CT))) %>%
#   dplyr::mutate(parasitesnum = as.numeric(as.character(parasites))) %>%
#   ggplot() +
#   geom_point(aes(x=parasitesnum, y=CTnum, shape=factor(parasitesnum), alpha=0.8), stat="identity") +
#   geom_hline(aes(yintercept=40), colour="#de2d26", linetype="dashed") +
#   ggtitle(label="Lower Limit of Detection for \n P. vivax Real-Time Assay") +
#   xlab("Parasites/uL") + ylab("Cycle Threshold Value") +
#   theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
#   theme(axis.title = element_text(hjust = 0.5, size=14)) +
#   theme(axis.text = element_text(hjust = 0.5, size=13)) +
#   theme(legend.position = "none")
# graphics.off()
