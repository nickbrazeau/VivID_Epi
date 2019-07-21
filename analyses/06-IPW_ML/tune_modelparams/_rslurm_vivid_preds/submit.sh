#!/bin/bash
#
#SBATCH --array=0-16
#SBATCH --job-name=vivid_preds
#SBATCH --output=slurm_%a.out

#SBATCH --mem=128000

#SBATCH --array=0-17%16

#SBATCH --cpus-per-task=8

#SBATCH --error=%A_%a.err

#SBATCH --output=%A_%a.out

#SBATCH --time=11-00:00:00
/nas/longleaf/apps/r/3.6.0/lib64/R/bin/Rscript --vanilla slurm_run.R
