#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="cellchat-env"

echo "==> Checking whether pkgbuild can find build tools (if R is already available)"
if command -v Rscript >/dev/null 2>&1; then
  Rscript -e "pkgbuild::check_build_tools(debug = TRUE)" || true
else
  echo "Rscript not found yet; skipping initial pkgbuild diagnostic"
fi

echo "==> Initialising conda/mamba shell support"
if ! command -v mamba >/dev/null 2>&1; then
  echo "Error: mamba is not in PATH"
  exit 1
fi

# Make 'mamba activate' work in non-interactive shells
eval "$(conda shell.bash hook)"

echo "==> Creating environment: ${ENV_NAME}"
mamba create -y -n "${ENV_NAME}" -c conda-forge r-base

echo "==> Activating environment: ${ENV_NAME}"
conda activate "${ENV_NAME}"

echo "==> Installing system/compiler and R package dependencies"
mamba install -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
  gcc_linux-64 \
  gxx_linux-64 \
  gfortran_linux-64 \
  make \
  pkg-config \
  cmake \
  autoconf \
  automake \
  libtool \
  binutils \
  zlib \
  libpng \
  cairo \
  pango \
  freetype \
  fontconfig \
  libjpeg-turbo \
  libxml2 \
  r-remotes \
  bioconductor-biocneighbors \
  r-igraph \
  r-sf \
  r-nmf \
  r-svglite \
  bioconductor-complexheatmap \
  bioconductor-biocgenerics \
  r-seurat \
  bioconductor-biobase \
  r-devtools \
  bioconductor-rhdf5

# echo "==> Installing SpatialCellChat GitHub dependencies"
# Rscript -e "remotes::install_github('JEFworks-Lab/MERINGUE', build_vignettes = TRUE)"
# Rscript -e "remotes::install_github('KlugerLab/ALRA')"
# Rscript -e "remotes::install_github('zdebruine/RcppML')"

# echo "==> Installing SpatialCellChat"
# Rscript -e "remotes::install_github('jinworks/SpatialCellChat', dependencies = TRUE)"

# echo "==> Testing SpatialCellChat installation"
# Rscript -e "library(SpatialCellChat)"

echo "==> Installing helper packages"
Rscript -e "remotes::install_github('cellgeni/schard')"
Rscript -e "remotes::install_github('immunogenomics/presto')"

echo "==> Done"


# new additions go here #

echo "==> Installing CellChat"
Rscript -e "remotes::install_github('jinworks/CellChat', dependencies = TRUE)"

echo "==> Testing CellChat installation"
Rscript -e "library(CellChat)"


# NOTE: Did not really end up using SpatialCellChat
# CellChat spatial functions work just fine
# There are clashes between the two packages
# and it is not worth the effort to resolve them at this time