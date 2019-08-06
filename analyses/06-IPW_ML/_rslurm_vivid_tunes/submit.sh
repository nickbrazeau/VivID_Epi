#!/bin/bash
#
#SBATCH --array=0-11
#SBATCH --job-name=vivid_tunes
#SBATCH --output=slurm_%a.out

#SBATCH --mem=32000

#SBATCH --array=0-12%13

#SBATCH --cpus-per-task=8

#SBATCH --error=%A_%a.err

#SBATCH --output=%A_%a.out

#SBATCH --time=11-00:00:00
/nas/longleaf/apps/r/3.6.0/lib64/R/bin/Rscript --vanilla slurm_run.R
