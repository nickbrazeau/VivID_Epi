library(Biostrings)
library(tidyverse)
devtools::install_github("mhahsler/rBLAST")

#-----------------------------
# import mtdna
#-----------------------------
pvp01 <- Biostrings::readDNAStringSet(
  filepath = "~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvP01.fasta"
)
pvp01 <- pvp01$PvP01_MIT_v1

pvsal1 <- Biostrings::readDNAStringSet(
  filepath = "~/Documents/MountPoints/mountIDEEL/resources/genomes/Pvivax/genomes/PvSal1.fasta"
)
pvsal1 <- pvsal1$`NC_007243 | organism=Plasmodium_vivax_Sal-1 | version=2013-05-01 | length=5990 | SO=mitochondrial_chromosome`




#-----------------------------
# import primers from PMID: 24639297
#-----------------------------
mtprm <- readxl::read_excel(path = "~/Documents/GitHub/VivID_Epi/WetLabWork/PvmtdDNADesign/pvmtdna_primers.xlsx", sheet = 1)

mtprm$sequencednastring <- lapply(mtprm$sequence, Biostrings::DNAString)

prm <- split(mtprm, f=factor(mtprm$primerorientation))

#-----------------------------
# match to pvsal1 and pvp01
#-----------------------------
primerlanding <- function(dnastring, target, orientation = NULL){
  # assumes 5' - 3' orientation
  # # returns 5'- 3' orientation
  if(orientation == "forward"){
    lnd <- Biostrings::matchPattern(dnastring, target)
    
    ret <- data.frame(sequence = as.character(lnd[[1]]),
                      start = lnd@ranges@start,
                      end = lnd@ranges@start + lnd@ranges@width - 1 # 1-based
    )
    
  } else if(orientation == "reverse"){
    lnd <- Biostrings::matchPattern(Biostrings::reverseComplement(dnastring), target)
    
    ret <- data.frame(sequence = as.character( Biostrings::reverseComplement(lnd[[1]]) ),
                      start = lnd@ranges@start,
                      end = lnd@ranges@start + lnd@ranges@width - 1 # 1-based
    )
  }
  # extract out info

  return(ret)
}

pvp01frwd <- lapply(prm$forward$sequencednastring, primerlanding, target = pvp01, orientation = "forward") %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(target = "pvp01") %>% 
  cbind.data.frame(step = prm$forward$step, .) # this ugliness because of the repeat forward primer

pvp01rvrs <- lapply(prm$reverse$sequencednastring, primerlanding, target = pvp01, orientation = "reverse") %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(target = "pvp01") %>% 
  cbind.data.frame(step = prm$reverse$step, .) # this ugliness because of the repeat forward primer

pvsal1frwd <- lapply(prm$forward$sequencednastring, primerlanding, target = pvp01, orientation = "forward") %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(target = "pvsal1") %>% 
  cbind.data.frame(step = prm$forward$step, .) # this ugliness because of the repeat forward primer


pvsal1rvrs <- lapply(prm$reverse$sequencednastring, primerlanding, target = pvp01, orientation = "reverse") %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(target = "pvsal1") %>% 
  cbind.data.frame(step = prm$reverse$step, .) # this ugliness because of the repeat forward primer


#-----------------------------
# visualize for PvP01
#-----------------------------
pvp01prm <- dplyr::bind_rows(pvp01frwd, pvp01rvrs) %>% 
  dplyr::left_join(x=mtprm, y = ., by = c("sequence", "step")) %>% 
  dplyr::select(c("primer", "step", "sequence", "primerorientation", "start", "end"))

pvp01prmplotdf <- pvp01prm %>% 
  dplyr::mutate(groupid = as.numeric(forcats::fct_rev(as.factor(primer))),
                groupid = (groupid)/20)  # ugly hack
pvp01prmplotdf <- pvp01prmplotdf %>% 
  dplyr::mutate(start_upd = ifelse(primerorientation == "reverse", end, start),
                end_upd = ifelse(primerorientation == "reverse", start, end),
                prmrorientstyle = stringr::str_to_title( paste(step, primerorientation) )
  )



ggplot() +
  geom_rect(aes(xmin = 0, xmax = length(pvp01)), ymin = 0, ymax = 0.8, fill = "#d9d9d9", color = "#d9d9d9") +
  geom_segment(data = pvp01prmplotdf,
               aes(x = start_upd, xend = end_upd, y = groupid, yend = groupid,
                   color = factor(prmrorientstyle)),
               arrow = arrow(length = unit(0.2,"cm"))) +
  scale_colour_manual("Primer \n Orientation", values = c("#a50026", "#4575b4", "#d73027", "#313695")) +
  ggtitle("Primer Distributions across PvP01 Genome for \n Outer and Inner Reactions") +
  xlab("P. vivax mtDNA Genome") + ylab("Primer Landing Sites") +
  theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 14), 
        axis.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 12), 
        axis.text = element_text(famil = "Arial", hjust = 0.5, size = 11), 
        legend.position = "right",
        legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.5, size = 12),
        legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank())
  

#-----------------------------
# visualize for PvSal1
#-----------------------------
pvsal1prm <- dplyr::bind_rows(pvsal1frwd, pvsal1rvrs) %>% 
  dplyr::left_join(x=mtprm, y = ., by = c("sequence", "step")) %>% 
  dplyr::select(c("primer", "step", "sequence", "primerorientation", "start", "end"))

pvsal1prmplotdf <- pvsal1prm %>% 
  dplyr::mutate(groupid = as.numeric(forcats::fct_rev(as.factor(primer))),
                groupid = (groupid)/20)  # ugly hack
pvsal1prmplotdf <- pvsal1prmplotdf %>% 
  dplyr::mutate(start_upd = ifelse(primerorientation == "reverse", end, start),
                end_upd = ifelse(primerorientation == "reverse", start, end),
                prmrorientstyle = stringr::str_to_title( paste(step, primerorientation) )
  )



ggplot() +
  geom_rect(aes(xmin = 0, xmax = length(pvsal1)), ymin = 0, ymax = 0.8, fill = "#d9d9d9", color = "#d9d9d9") +
  geom_segment(data = pvsal1prmplotdf,
               aes(x = start_upd, xend = end_upd, y = groupid, yend = groupid,
                   color = factor(prmrorientstyle)),
               arrow = arrow(length = unit(0.2,"cm"))) +
  scale_colour_manual("Primer \n Orientation", values = c("#a50026", "#4575b4", "#d73027", "#313695")) +
  ggtitle("Primer Distributions across PvSal1 Genome for \n Outer and Inner Reactions") +
  xlab("P. vivax mtDNA Genome") + ylab("Primer Landing Sites") +
  theme(plot.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 14), 
        axis.title = element_text(famil = "Arial", face = "bold", hjust = 0.5, size = 12), 
        axis.text = element_text(famil = "Arial", hjust = 0.5, size = 11), 
        legend.position = "right",
        legend.title = element_text(famil = "Arial", face = "bold", vjust = 0.5, size = 12),
        legend.text = element_text(famil = "Arial", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank())



