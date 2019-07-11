

# Going to use the primer sets from PMC5789620
# I adapted the Forward Primer not overlap with the
# Ape Phylogenetic SNP
# My Forward: 5’-TTATATCCACCATTAAGTACATCACTT-3’
# Original Forward: 5’-AACCTTTAGATCTTAGATGCATTACA-3’
# Reverse: 5’-AGTGTTAAACCTTTAGATCTTAGATG-3’


# Now need to design unique MIDs that are far away from the primer set and 
# follow all of the amplicon rules



# imports
library(tidyverse)
library(Biostrings)
library(seqinr)
remotes::install_github("IDEELResearch/SeekDeepRANN")
library(SeekDeepRANN)

#-----------------------------
# import Pv mtdna
#-----------------------------
pvp01 <- Biostrings::readDNAStringSet("~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta", "fasta")
pvmtdna <-  pvp01[ grepl("PvP01_MIT_v1", pvp01@ranges@NAMES) ]

# forward landing
Biostrings::vmatchPattern("TTATATCCACCATTAAGTACATCACTT", pvmtdna) 
# reverse landing
Biostrings::vmatchPattern(Biostrings::reverseComplement(Biostrings::DNAString("AGTGTTAAACCTTTAGATCTTAGATG")), pvmtdna)
# target gene
seqinr::c2s( seqinr::s2c(Biostrings::toString(pvmtdna))[3726:3880] )


# now read in design fasta for the purpose of finding MIDs
# and find MIDs

mids <- SeekDeepRANN::MIDPrimerFinder(
  design_fasta = "WetLabWork/Pv_MIDs/mid_design_fasta.fa",
  MID2MIDmatchesAllowed = 4,
  MIDlength = 8,
  MIDhomopolymerallowance = 2,
  MID2targetmismatchesAllowed = 3,
  MIDnum = 20
)


saveRDS(mids, "WetLabWork/Pv_MIDs/final_mids.rds")
