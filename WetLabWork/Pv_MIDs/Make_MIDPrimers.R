

library(tidyverse)
library(Biostrings)


designfa <- Biostrings::readDNAStringSet("WetLabWork/Pv_MIDs/mid_design_fasta.fa")

forwardprimer <- designfa$` forward`
reverseprimer <- designfa$` reverse`


# confirm one last time that everything is in order
Biostrings::matchPattern(forwardprimer, designfa$` cox1`)
Biostrings::matchPattern(Biostrings::reverseComplement(reverseprimer), designfa$` cox1`)

# our amplicon including our Primers is 155-bp. With MIDs of 8 bp on both sides, we have a 171 bp fragment to sequence 

# Pull in MIDs
mids <- readRDS(file = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/WetLabWork/Pv_MIDs/final_mids.rds")

# to order
fw <- paste0(mids$finalMIDs[1:10], as.character(Biostrings::BString(forwardprimer)))
rv <- paste0(mids$finalMIDs[11:20], as.character(Biostrings::BString(reverseprimer)))
to.order <- data.frame( 
  name = c(paste0("PvcoxI_F_M", 1:10), paste0("PvcoxI_R_M", 11:20)),
  primer = c(fw, rv)
    )
write_csv(x = to.order, path = "~/Desktop/Pvmtcox-Primers.csv")


fw <- paste0("fw", 1:12)
rv <- paste0("rv", 13:20)
grid <- expand.grid(fw,rv)
write_csv(x = grid, path = "~/Desktop/primermap.csv")

