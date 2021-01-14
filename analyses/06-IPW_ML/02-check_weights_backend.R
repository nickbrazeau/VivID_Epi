#----------------------------------------------------------------------------------------------------
# Purpose of this script is to check for 
# covariate balance after applying the weights
#----------------------------------------------------------------------------------------------------
library(drake)
library(tidyverse)
library(energy)
source("R/00-IPTW_functions.R")
set.seed(48, "L'Ecuyer")

#............................................................
# Import Data and superlearner results
#............................................................
dt <- readRDS("data/derived_data/vividepi_recode_completecases.rds")
params <- readRDS("results/ensembl_cvRisk_paramdf.RDS")

#............................................................
# Make Model Map
#............................................................
select_samples_by_iptw_and_feattarget <- function(target, feat, iptw, data){
  wi <- sample(1:nrow(data), prob = iptw, replace = T)
  ret <- data %>% 
    dplyr::select(c(target, feat))
  
  ret <- ret[wi, ]
  return(ret)
  
}


#............................................................
# Note, I want to sample 100 realizations of my
# sample weights since I can't apply the weights directly
# in the energy function 
#............................................................
modelmap.ml <- params %>% 
  dplyr::select(c("task", "iptw", "target"))

modelmap.ml <- replicate(n = 100, modelmap.ml,
                         simplify = F) %>% 
  dplyr::bind_rows()

modelmap.ml$feat <- map(modelmap.ml$task, mlr::getTaskFeatureNames)
modelmap.ml$data <- lapply(1:nrow(modelmap.ml), function(x) return(dt))

# unnest to get all feats
modelmap.ml <- modelmap.ml %>% 
  tidyr::unnest(cols = feat)



modelmap.ml$data <- purrr::pmap(modelmap.ml[,c("target", "feat", "iptw", "data")], 
                                select_samples_by_iptw_and_feattarget)

modelmap.ml <- modelmap.ml %>% 
  dplyr::select(c("target", "feat", "data")) %>% 
  dplyr::rename(covar1 = target,
                covar2 = feat) %>% 
  dplyr::mutate(lvl = "IPTW")


# get base target-feats too
modelmap.nowi <- params %>% 
  dplyr::select(c("task", "target")) 

modelmap.nowi$feat <- map(modelmap.nowi$task, mlr::getTaskFeatureNames)
modelmap.nowi$data <- lapply(1:nrow(modelmap.nowi), function(x) return(dt))

modelmap.nowi <- modelmap.nowi %>% 
  tidyr::unnest(col = feat) %>% 
  dplyr::select(c("target", "feat", "data")) %>% 
  dplyr::rename(covar1 = target,
                covar2 = feat) %>% 
  dplyr::mutate(lvl = "base")



modelmap <- rbind.data.frame(modelmap.nowi, modelmap.ml)


#............................................................
# Energy function for slurm
#............................................................

slurm_calc_corr <- function(covar1, covar2, data, lvl, id){
  x1 <- unlist( data[,covar1] )
  x2 <- unlist( data[,covar2] )
  if(is.factor(x1)){
    x1 <- as.numeric(x1)
  }
  if(is.factor(x2)){
    x2 <- as.numeric(x2)
  }
  
  ret <- energy::dcor(x = x1, y = x2)
  
  return(ret)
}


#............................................................
# parallelize with drake
#...........................................................
workers <- 512 # nodes to ask for, fewer nodes, less expensive for reading in data and not placing burden on scheduler
library(drake)

#......................
# functions for drake
#......................
drake_wrapper <- function(batchset_df) {
  
  # call future
  no_cores <- future::availableCores() - 1
  if (no_cores > 1) {
    future::plan(future::multicore, workers = no_cores)
  } else {
    future::plan("sequential")
  }
  
  
  #......................
  # run batches to not overload scheduler
  #......................
  batchset_df$energcorr <- furrr::future_pmap(batchset_df, slurm_calc_corr)
  
  # now write out
  dir.create("/pine/scr/n/f/nfb/Projects/VivID_Epi/06-IPW_ML/wi_covar_corr/", recursive = T)
  saveRDS(batchset_df,
          file = paste0("/pine/scr/n/f/nfb/Projects/VivID_Epi/06-IPW_ML/wi_covar_corr/",
                        "batchset_", unique(batchset_df$id), ".RDS"))
  
}

#......................
# make batches
#......................
batchnum <- sort( rep(1:workers, ceiling(nrow(modelmap) / workers)) )
batchnum <- batchnum[1 :nrow(modelmap)]

modelmap_nested <- modelmap %>%
  dplyr::mutate(id = batchnum,
                batchset = batchnum) %>%
  dplyr::group_by(batchset) %>%
  tidyr::nest() %>%
  dplyr::ungroup()


#......................
# make drake plan
#......................
batch_names <- paste0("batch", modelmap_nested$batchset)
plan <- drake::drake_plan(
  runs = target(
    drake_wrapper(data),
    transform = map(
      .data = !!modelmap_nested,
      .names = !!batch_names
    )
  )
)

#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_clst/slurm_clustermq_LL_short.tmpl")
make(plan,
     parallelism = "clustermq",
     jobs = nrow(modelmap_nested),
     log_make = "wi_covar_ipw_drake.log", verbose = 4,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     garbage_collection = TRUE,
     lock_envir = FALSE,
     lock_cache = FALSE)
