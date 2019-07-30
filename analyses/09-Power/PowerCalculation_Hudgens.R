library(pwr)


## NOTES


powercalculator.glmRR <- function(n=15879, exp_prob=0.5, p=0.05, p0=0.02){
  
  df <- data.frame(obs=factor(seq(1:n)),
                   exp=sample(x=c(0,1), size=n, replace = T, prob=c(exp_prob, 1-exp_prob))) # df of exposure
  p <- 2*p # inv average prev for both groups 
  p0 <- p0 # prev among unexposed
  p1 <- p-p0 # prev among exposed
  RR <- exp(log(p1) - log(p0))


    df$dz[df$exp == 1] <- rbinom(sum(df$exp == 1),1,p1)
    df$dz[df$exp == 0] <- rbinom(sum(df$exp == 0),1,p0)

    mod <- glm(dz ~ exp, data=df,
                  family=binomial(link="log"))

    pi <- broom::tidy(mod)$p.value[2]

    ret <- data.frame(RR=RR, p=pi)

    return(ret)

}



### run lots of these at different levels of p0
p0sim <- seq(0.01, 0.032, by=0.001)
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
