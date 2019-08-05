#!/bin/bash

<<<<<<< HEAD
#SBATCH --ntasks=16
=======
#SBATCH --ntasks=1
>>>>>>> master
#SBATCH --time=5-00:00:00
#SBATCH --mem=49512
#SBATCH --mail-type=all
#SBATCH --mail-user=nbrazeau@med.unc.edu

<<<<<<< HEAD
Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi"); source("analyses/09-Power/01-PowerCalculation_Hudgens.R")'
=======
Rscript -e 'setwd("/proj/ideel/meshnick/users/NickB/Projects/VivID_Epi"); source("analyses/09-Power/01-_PowerCalculation_Hudgens.R")'
>>>>>>> master
