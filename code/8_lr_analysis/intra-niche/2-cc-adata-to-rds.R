# the purpose of this script is to convert the h5ad file created in the previous 
# step into an RDS file that can be easily loaded into R for downstream analysis.

library(Seurat)
library(schard)

MAIN_DIR <- "/mnt/data/project0062/proseg_data"
setwd(MAIN_DIR)

## Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

## Object versions / paths
CUR_OBJ_V <- "annotated"
CC_OBJ_PATH  <- file.path(
  MAIN_DIR, "data", "comb-full", "h5ad", sprintf("comb-full-%s-cellchat.h5ad", CUR_OBJ_V)
)

# load h5ad as Seurat
comb = schard::h5ad2seurat(CC_OBJ_PATH)

# Remove cells labelled "Amb" in ann_lvl_3_refined
comb <- subset(
  comb,
  subset = ann_lvl_3_refined != "Amb"
)

# NOTE: The counts and data layers contain the same values

# %% Save the comb object as RDS 

CC_RDS_PATH <- file.path(
  MAIN_DIR, "data", "comb-full", "rds", sprintf("comb-full-%s-cellchat.rds", CUR_OBJ_V)
)

# create directory if it doesn't exist
if (!dir.exists(dirname(CC_RDS_PATH))) {
  dir.create(dirname(CC_RDS_PATH), recursive = TRUE)
}

saveRDS(comb, CC_RDS_PATH)
