library(Biostrings)
library(tidyverse)
remotes::install_github("mhahsler/rBLAST")
library(rBLAST)


#-----------------------------
# import Pv mtdna
#-----------------------------

pvp01 <- Biostrings::readDNAStringSet("~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta", "fasta")
pvp01 <-  pvp01[ grepl("PvP01_MIT_v1", pvp01@ranges@NAMES) ]
vmatchPattern("TTATATCCACCATTAAGTACATCACTT", pvp01) 
vmatchPattern(reverseComplement(DNAString("AGTGTTAAACCTTTAGATCTTAGATG")), pvp01)



#-----------------------------
# import other mtdna
#-----------------------------
# pf3d7
pf3d7 <- Biostrings::readDNAStringSet("~/Documents/MountPoints/mountIDEEL/resources/genomes/Pfalciparum/genomes/Pf3D7.fasta", "fasta")
pf3d7 <- pf3d7[ grepl("M76611", pf3d7@ranges@NAMES) ]

vmatchPattern("CAGCAGCAGAATTTGGAGGAG", pf3d7, max.mismatch = 2)
vmatchPattern(reverseComplement(DNAString("ATGAGCCCATACAACACTTCC")), pf3d7, max.mismatch = 3)




malariae <- Biostrings::readDNAStringSet("WetLabWork/PvAmplicon/PvmtdDNADesign/data/PlasmoDB-38_PmalariaeUG01_Genome.fasta", "fasta")
malariae <- malariae[grepl("PmUG01_MIT_v1", malariae@ranges@NAMES) ]

vmatchPattern("CAGCAGCAGAATTTGGAGGAG", malariae, max.mismatch = 2)
vmatchPattern(reverseComplement(DNAString("ATGAGCCCATACAACACTTCC")), malariae, max.mismatch = 2)

ovale <- Biostrings::readDNAStringSet("WetLabWork/PvAmplicon/PvmtdDNADesign/data/PlasmoDB-38_PovalecurtisiGH01_Genome.fasta", "fasta")
#ovale <- ovale[grepl("PmUG01_MIT_v1", ovale@ranges@NAMES) ]

vmatchPattern("CAGCAGCAGAATTTGGAGGAG", ovale, max.mismatch = 2)
vmatchPattern(reverseComplement(DNAString("ATGAGCCCATACAACACTTCC")), ovale, max.mismatch = 2)




# GenBank: MF197850.1
# Plasmodium simium isolate
psim <- Biostrings::readDNAStringSet("WetLabWork/PvAmplicon/PvmtdDNADesign/data/MF197850_1.fasta")

# GenBank: KX645965.1
# Plasmodium cf. inui
pinfui <- Biostrings::readDNAStringSet("WetLabWork/PvAmplicon/PvmtdDNADesign/data/KX645965_1.fasta")

# GenBank: KY790474.1
# Plasmodium gaboni
pinfui <- Biostrings::readDNAStringSet("WetLabWork/PvAmplicon/PvmtdDNADesign/data/KY790474_1.fasta")






# Can't tell the difference to psim but its a NMW isolate... not OWM?
toString(subseq(mtdna, 3726, 3880))
toString(subseq(psim, 3727, 3881))


toString(subseq(pinfui, 2985, 3880))















# --------------------------------
# Previous code chunks that may be useful
library(rentrez)






#####################################################
#################  Locate Data    ##################
#####################################################
# Am going to query the ncbi database for vivax mtDNA 

entrez_db_searchable("nucleotide")

#####################################
#####    Download mtDNA    #########
####################################
vivaxmtdna.mtDNA.search <- entrez_search(db="nucleotide", 
                                         term="((plasmodium vivax[ORGN] OR Plasmodium vivax[ORGN]  AND (mitochondrion OR mitochondrial OR mitochondria OR cytb OR coxI OR coxIII))",
                                         retmax=10000,
                                         use_history=T)

vivaxmtdna.mtDNA.search


# initialize metadatadb
vivaxmtdna.mtDNA.mtdt <- as.data.frame(matrix(ncol=6, nrow=length(vivaxmtdna.mtDNA.search$ids)))
colnames(vivaxmtdna.mtDNA.mtdt) <- c("organism", "title", "accessionnum", "subnamecountry", "completeness", "seqencelength")  


#for(seq.start in seq(1,length(vivaxmtdna.mtDNA.search$ids),1)){
for(sq in seq(1,length(vivaxmtdna.mtDNA.search$ids),1)){
  # Pulling down the fasta sequences 
  fastas <- entrez_fetch(db="nucleotide", rettype="fasta",
                         id = vivaxmtdna.mtDNA.search$ids[sq])
  cat(fastas, file=paste0("~/Documents/MountPoints/mountedMeshnick/Projects/DRC_Vivax/World_Pv_Sequences/mtDNA_Analysis/", Sys.Date(), "pvmtdnaNucleotidesearch.fasta"), append=TRUE)
  cat(sq+1, "sequences downloaded\r")
  # Pulling down the relevant metadata
  mtdt <- entrez_summary(db="nucleotide", id = vivaxmtdna.mtDNA.search$ids[sq])
  
  
  vivaxmtdna.mtDNA.mtdt[sq,1] <-  mtdt$organism
  vivaxmtdna.mtDNA.mtdt[sq,2] <-  mtdt$title
  vivaxmtdna.mtDNA.mtdt[sq,3] <-  mtdt$accessionversion
  vivaxmtdna.mtDNA.mtdt[sq,4] <-  mtdt$subname
  vivaxmtdna.mtDNA.mtdt[sq,5] <- mtdt$completeness
  vivaxmtdna.mtDNA.mtdt[sq,6] <- mtdt$slen
  
}


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4047736/