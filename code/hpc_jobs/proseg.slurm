#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=SegCPU
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1   # one process per GPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-user=robin.carvajal@glasgow.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# Load cargo environment (same as your .bashrc)
. "/mnt/data/project0062/.cargo/env"

# Optional: confirm Rust is active
echo "Cargo version: $(cargo --version)"
echo "Rustc version: $(rustc --version)"

# === Run ProSeg ===

# proseg \
#   --cosmx /mnt/data/project0062/cosmx_gray/data/proseg_samples/CTRL01/CTRL01_tx_file.csv.gz \
#   --nthreads 64 \
#   --voxel-layers 1 \
#   --voxel-size 0.5 --burnin-voxel-size 1.0 \
#   --cell-compactness 0.03 \
#   --diffusion-probability 0.02 \
#   --nuclear-reassignment-prob 0.05 \
#   --max-transcript-nucleus-distance 60 \
#   --samples 500 --burnin-samples 1000 \
#   --output-spatialdata /mnt/data/project0062/cosmx_gray/data/proseg_samples/CTRL01/output.zarr \
#   --overwrite

proseg \
  --cosmx /mnt/data/project0062/cosmx_gray/data/proseg_samples/CTRL01/CTRL01_tx_file.csv.gz \
  --nthreads 64 \
  --voxel-layers 1 \
  --voxel-size 1.0 --burnin-voxel-size 2.0 \
  --cell-compactness 0.03 \
  --diffusion-probability 0.02 \
  --nuclear-reassignment-prob 0.01 \
  --max-transcript-nucleus-distance 60 \
  --samples 500 --burnin-samples 1000 \
  --output-spatialdata /mnt/data/project0062/cosmx_gray/data/proseg_samples/CTRL01/output.zarr \
  --overwrite
