#!/bin/bash
#
#SBATCH --array=0-9
#SBATCH --job-name=vivid_spSL_final
#SBATCH --output=slurm_%a.out

#SBATCH --mem=32000

#SBATCH --array=0-10%17

#SBATCH --cpus-per-task=1

#SBATCH --error=%A_%a.err

#SBATCH --output=%A_%a.out

#SBATCH --time=5-00:00:00
/nas/longleaf/apps/r/3.6.0/lib64/R/bin/Rscript --vanilla slurm_run.R
