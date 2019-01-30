#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to visualize potential primer sites
# for Pvamplicon work
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(vcfR)

#.................................
# read in
#.................................
# note, I subsetted this vcf on the command line previously
dv <- vcfR::read.vcfR(file = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/WetLabWork/PvampliconDesign/vcf/prop-primers_snps.VQSR900.human_only.vcf.gz")
bd <- readr::read_tsv(file = "WetLabWork/primerbeds/pv_p01_prop-primers.bed")

#.................................
# datawrangle
#.................................
bd <- bd %>%
  dplyr::rename(seqname = chr)
bdlist <- split(bd, factor(1:nrow(bd)))

bd$vcfRobj <- purrr::map(bdlist, vcfRmanip::vcfR2SubsetChromPos, vcfRobject = dv)
bd$nvar <- unlist(purrr::map(bd$vcfRobj, function(x){ return(nrow(x@gt)) }))

#.................................
# data manip for plotting
#.................................

bd$plotdf <- purrr::map(bd$vcfRobj, function(x, metadata = mt){
                              ret <- cbind.data.frame(CHROM = vcfR::getCHROM(x),
                                     POS = vcfR::getPOS(x),
                                     vcfR::extract.gt(x, element = "GT")) %>% 
                                     gather(., key = "iid", value = "gt", 3:ncol(.)) %>%
                                     dplyr::mutate(gt_fct = factor( gt, levels = c("0/0", "0/1", "1/1"), 
                                                                    labels = c("REF", "HET", "ALT")),
                                                                    POS = as.numeric(POS)) %>% 
                                     inner_join(metadata, ., by = "iid") %>% 
                                     dplyr::mutate(cont = rplasmodium::encode_continent(country)) %>% 
                                     dplyr::mutate(country_fct = factor(country)) 
                              
                              return(ret)
              
})


#.................................
# Make Plots
#.................................
bd$vcfplots <- purrr::map2(bd$plotdf, bd$geneid, function(x, y){
                          ret <-  x %>% 
                                  dplyr::arrange(POS) %>% 
                                  ggplot() +
                                  geom_tile(aes(x=factor(iid), y= factor(POS), fill = gt_fct)) + 
                                  facet_wrap(cont + country_fct ~ ., scales = "free_x", nrow = 1) +
                                  ggtitle(paste("Nucleotide variation for Gene", y)) +
                                  theme_classic() +
                                  theme(axis.text.x = element_blank(),
                                        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
                          return(ret)
  
})














