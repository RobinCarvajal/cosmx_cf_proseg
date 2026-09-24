# The purpose of this script is to create separate Seurat objects for each niche 
# defined in the metadata of the combined object.


# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


# Set paths --------------------------------------------------------------------
MAIN_DIR <- "/mnt/data/project0062/proseg_data"
setwd(MAIN_DIR)

# importing cellchat functions
source("code/8_lr_analysis/cellchat-functions.R")

# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Object versions / paths
CUR_OBJ_V <- "annotated"
CC_OBJ_PATH  <- file.path(
  MAIN_DIR, "data", "comb-full", "rds", sprintf("comb-full-%s-cellchat.rds", CUR_OBJ_V)
)

# where to save the objects
OUT_DIR <- file.path(
  MAIN_DIR, "results", "comb-full", "lr-manual-niches"
)

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# Load RDS as Seurat -----------------------------------------------------------
comb = readRDS(CC_OBJ_PATH)

# Separate the objet in subsets by niches and then condition -------------------

NICHES_COL <- "niche_name"

meta_df <- comb@meta.data
niche_names <- unique(meta_df[[NICHES_COL]])

for (niche in niche_names) {

  message(sprintf("Working on: %s", niche))

  filter_mask <- meta_df[[NICHES_COL]] == niche
  cells_use <- rownames(meta_df)[filter_mask]

  niche_obj <- subset(comb, cells = cells_use)

  CUR_NICHE_DIR <- file.path(OUT_DIR, as.character(niche))
  if (!dir.exists(CUR_NICHE_DIR)) {
    dir.create(CUR_NICHE_DIR, recursive = TRUE)
  }

  # save RDS 
  saveRDS(niche_obj, file.path(CUR_NICHE_DIR, "niche-obj.rds"))

}