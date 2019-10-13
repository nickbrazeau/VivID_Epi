#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=16g
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi"); source("analyses/07-spatial_prediction/03-Pv_get_predictions_backend.R")'