#----------------------------------------------------------------------------------------------------
# Purpose of this script is analyze the Duffy Sanger Sequencing Data
# generated on the Pv cases
#----------------------------------------------------------------------------------------------------
# libraries
library(tidyverse)
library(Biostrings)

# functions
readseq <- function(filepath){
  
  nm <- tolower(base::basename(filepath))
  nm <- gsub(pattern = "vivax", replacement = "", x = nm)
  nm <- gsub(pattern = "prime", replacement = "", x = nm)
  ret <- data.frame(
    smpl = stringr::str_match(nm, "[:alnum:][:alnum:][:alnum:][:alnum:][:alnum:]"), # insane naming scheme from sequencer, thankfully only sample barcode has 5 alphanumeric
    ornt = ifelse(stringr::str_detect(nm, "for|fw"), "forward", "reverse"),
    seq = sapply(read.table(filepath, header=F), paste, collapse = "")
    )
  
  ret$smpl <- ifelse(is.na(ret$smpl), stringr::str_match(nm, "pc1"), as.character(ret$smpl))
  ret$smpl <- ifelse(is.na(ret$smpl), stringr::str_match(nm, "pc2"), as.character(ret$smpl)) 
  
  return(ret)
}

seq2fasta <- function(seqdf, outpath){
  cat(file=paste0(outpath, "/", seqdf$smpl, "-", seqdf$ornt, ".fa"), 
      paste(paste0(">", seqdf$smpl, "-", seqdf$ornt),
            seqdf$seq, sep="\n")
  )
  return(0)
}


#---------------------------------------------------
# Read in seq and convert to fasta
#---------------------------------------------------
drs <- list.dirs(path = "/Volumes/share/1. Data/1. Raw Data/Adult_PvDARC_Genotyping/Sequences/", 
          full.names = T)[2:10]

fls <- unlist(lapply(drs, function(x){
  list.files(path = x, pattern = ".seq", full.names = T)
  }))

seq <- lapply(fls, readseq) %>% dplyr::bind_rows(.)


# convert the .seq files into fasta files
split(seq, 1:nrow(seq)) %>% 
  lapply(., seq2fasta, outpath = "/Volumes/share/1. Data/1. Raw Data/Adult_PvDARC_Genotyping/fastas_from_R/")


# Primers are from [Menard et al. 2010](http://www.pnas.org/content/pnas/suppl/2010/03/12/0912496107.DCSupplemental/pnas.200912496SI.pdf). **These are the inner primers from a hemi-nested PCR protocol**. 

# note these are nested primers
Fw <- Biostrings::DNAString("GTGGGGTAAGGCTTCCTGAT") # 5'-3'
Rv <-  Biostrings::reverseComplement(Biostrings::DNAString("CAAACAGCAGGGGAAATGAG"))

target <- Biostrings::readDNAStringSet(file = "~/Documents/GitHub/VivID_Epi/hostgenotyping/DARC_analyses/X85785_1.fasta", 
                                       format="fasta", use.names = T)
lapply(list(Fw, Rv), function(x) {
  vmatchPattern(x, target)
})

target <- target[[1]][367:589] # from where primers land above

hostgnfor <- lapply(list.files("/Volumes/share/1. Data/1. Raw Data/Adult_PvDARC_Genotyping/fastas_from_R/",
                    full.names = T, pattern = "forward"), 
                function(x){Biostrings::readDNAStringSet(x, format = "fasta")})

hostgnrev <- lapply(list.files("/Volumes/share/1. Data/1. Raw Data/Adult_PvDARC_Genotyping/fastas_from_R/",
                               full.names = T, pattern = "reverse"), 
                    function(x){Biostrings::readDNAStringSet(x, format = "fasta")})



# cat(file = "/Volumes/share/1. Data/1. Raw Data/Adult_PvDARC_Genotyping/fastas_from_R/DARC_host_genomes_combined-forward.fa",
#     sapply(hostgnfor, names)
# 
# 
# names(hostgnfor[[1]])
# as.character(hostgnfor[[1]])
# 
# cat(file=paste0(getwd(), "/data/ama1hapsfaseg.fa"), paste(paste0(">", names(ama1hapsfaseg)),
#                                                           sapply(ama1hapsfaseg, paste, collapse=""), sep="\n"), sep="\n")
# t <- pairwiseAlignment(target, hostgnfor[[1]])
