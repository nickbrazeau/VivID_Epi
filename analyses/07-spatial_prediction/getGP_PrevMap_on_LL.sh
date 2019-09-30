#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi"); source("analyses/07-spatial_prediction/02-Pv_clst_map_GaussianProc_MCMC_backend.R")'