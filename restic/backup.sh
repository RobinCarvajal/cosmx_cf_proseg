#!/usr/bin/env bash
#SBATCH --account=project0062
#SBATCH --job-name=ResticBackup
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/data/project0062/cosmx_cf_proseg/restic/logs/%x-%j.out
#SBATCH --error=/mnt/data/project0062/cosmx_cf_proseg/restic/logs/%x-%j.err
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load apps/miniforge
conda activate restic-env

project_dir="/mnt/data/project0062/cosmx_cf_proseg" # absolute path is better
tag_name="$(basename "$project_dir")"
host_name="mars"

restic backup "$project_dir" --tag "$tag_name" --host "$host_name"