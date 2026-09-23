#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account="project0062"
#SBATCH --job-name=MergeAndConvert
#SBATCH --partition=nodes
#SBATCH --time=0-08:00:00
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=ALL

set -euo pipefail

############# LOADING MODULES #############
module load apps/miniforge

############ DIRS ###############
scripts_dir="$P62/cosmx_gray/scripts/1_pre_analysis/"

############ ACTIVATE CONDA ENV ###############

conda activate r-cosmx

############ CODE ##############

# MERGE
#Rscript "${scripts_dir}/1_merging_slides.R"

# CONVERT
Rscript "${scripts_dir}/2_rds_to_h5ad.R"

