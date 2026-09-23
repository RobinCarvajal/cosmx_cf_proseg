# %% Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(pryr) # to check ram usage mem_used()

# %% Paths

data_dir <- "/mnt/data/project0062/cosmx_gray/data"

rds_path <- file.path(data_dir, "combined/rds/Combined_v0.RDS")
h5seurat_path <- file.path(data_dir, "combined/h5ad/Combined_v0.h5Seurat")

# %% read obj.rds
obj_rds <- readRDS(rds_path)

# %% Keep only necessary data in the object

# Clear the images slot
obj_rds@images <- list()

# Diet the object
diet_obj <- DietSeurat(
  obj_rds,
  assays = "RNA",
  layers = "counts"
)

# export as "h5Seurat" object in your current working directory
SeuratDisk::SaveH5Seurat(
  diet_obj,
  filename = h5seurat_path,
  overwrite = TRUE
)

# convert to h5ad format
SeuratDisk::Convert(
  h5seurat_path,
  dest = "h5ad",
  overwrite = TRUE
)

# remove the temporary .h5seurat file
file.remove(h5seurat_path)
