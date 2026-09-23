#!/bin/bash
#SBATCH --account=project0062
#SBATCH --job-name=SegCPU
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1   # one process per GPU
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=48:00:00
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


SLIDE_NAMES=(SLIDE01 SLIDE02 SLIDE03) # SLIDE04 has been done

BASE_DIR="/mnt/data/project0062/reseg_test"

for slide in "${SLIDE_NAMES[@]}"; do
  tx="${BASE_DIR}/${slide}_tx_file.csv.gz"
  out="${BASE_DIR}/${slide}_proseg.zarr"

  echo "==> ${slide}"
  echo "  input : ${tx}"
  echo "  output: ${out}"

  # fail early if input missing
  [[ -r "$tx" ]] || { echo "ERROR: missing/unreadable $tx" >&2; exit 1; }

  proseg --cosmx "$tx" --output-spatialdata "$out" 
done