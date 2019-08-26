#--------------------------------------------------------------------------
# Purpose of this script is to collect the null distribution
# in order to use it for our performance measure
#--------------------------------------------------------------------------
# dependencies
library(tidyverse)

#...................................
# Import Data
#...................................
null.ret.files <- list.files(path = "~/Documents/MountPoints/mountedMeshnick/Projects/VivID_Epi/analyses/06-IPW_ML/00-null_distributions/_rslurm_nulldist_vivid_covar/",
                             pattern = ".RDS", full.names = T)

# I don't need the function or param file here
null.ret.files <- null.ret.files[ !basename(null.ret.files) %in% c("params.RDS", "f.RDS") ]
null.ret <- lapply( null.ret.files, function(x){ 
                                                ret <- readRDS(x)
                                                ret <- data.frame(
                                                  covar = names(ret),
                                                  corr = unlist(ret)
                                                )
                                                } ) %>% 
  dplyr::bind_rows()



null.ret <- null.ret %>% 
  group_by(covar) %>% 
  tidyr::nest() %>% 
  dplyr::rename(target = covar,
                nulldist = data)


saveRDS(object = null.ret, file = "analyses/06-IPW_ML/00-null_distributions/null_dist_return.RDS")
