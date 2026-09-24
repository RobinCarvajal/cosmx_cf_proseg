# The purpose of this script is to create CellChat objects for each condition 
# in the dataset, which will be used for ligand-receptor interaction analysis.

# this is the "cli" version of the script, which is meant to be run on the command line

# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
#library(SpatialCellChat)
library(CellChat)
library(patchwork)
library(future)
options(stringsAsFactors = FALSE)
library(argparse)

# Set up argument parser -------------------------------------------------------
parser <- ArgumentParser(
  description = "Create CellChat objects for each condition in the dataset."
)

parser$add_argument(
  "--main_dir", 
  type = "character", 
  default = "/data/cosmx_gray", 
  help = "Main directory path"

)
parser$add_argument(
  "--input_obj_path", 
  type = "character", 
  default = "results/comb/lr-manual-niches/1/niche-obj.rds", 
  help = "Path to the input niche object RDS file"
)

parser$add_argument(
  "--output_dir", 
  type = "character", 
  default = "results/comb/lr-manual-niches/1", 
  help = "Directory to save the output CellChat objects"
)

parser$add_argument(
  "--cc_functions_path", 
  type = "character", 
  default = "code/8_lr_analysis/cellchat-functions.R", 
  help = "Path to the CellChat functions script"
)

args <- parser$parse_args()


# Set paths --------------------------------------------------------------------
MAIN_DIR <- args$main_dir
setwd(MAIN_DIR)

cellchat_functions_path <- args$cc_functions_path
if (!file.exists(cellchat_functions_path)) {
  stop(
    paste(
      "The cellchat functions script was not found at the specified path:", 
      cellchat_functions_path
    )
  )
}
source(cellchat_functions_path)

# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Set up parallel processing ---------------------------------------------------
future::plan("sequential")

options(
  future.globals.maxSize = 120 * 1024^3,
  future.rng.onMisuse = "ignore"
)
# helps speed up LR identification

# Niche object path ------------------------------------------------------------

IN_OBJ_PATH <- file.path(MAIN_DIR, args$input_obj_path)

OUT_DIR <- file.path(MAIN_DIR, args$output_dir)

# load the niche object
niche_obj <- readRDS(IN_OBJ_PATH)

# Considerations for CosMx -----------------------------------------------------

# This conversion factor is used to take CosMx from pixels to micrometers
CONVERSION_FACTOR <- 0.12028
# This contact rannge is given in micrometers, it's standard for human cells
CONTACT_RANGE <- 10
# interaction range = the maximum interaction/diffusion length of ligands, in microns
INTERACTION_RANGE <- 250 # default 
# scale distance = the power of distance used to scale the interaction probability
# default is 0.01 but sometimes it ask for higher values
# adjust as neccessary
SCALE_DISTANCE <- 1

# Defining the database --------------------------------------------------------
cellchatDB <- CellChatDB.human
# keep only protein interactions
cellchatDB.use <- subsetDB(cellchatDB) 
# default is all the database, except non-protein interactions


# Separate the objet in subsets by niches and then condition -------------------

CONDITION_COL <- "condition"
CELLTYPE_COL <- "ann_lvl_3_refined"
NICHE_COL <- "cellcharter_niches"
SAMPLE_COL <- "sample_name"
GROUP_COLS <- c(CONDITION_COL, NICHE_COL, CELLTYPE_COL)
COORD_COLS <- c("centroid_x", "centroid_y")

# Get condition names without splitting the full object
cond_names <- unique(niche_obj@meta.data[[CONDITION_COL]])

# Process one condition at a time ----------------------------------------------

for (cond in cond_names) {

  message("Processing condition: ", cond)

  # Subset only the current condition
  cond_obj <- subset(
    niche_obj,
    subset = condition == cond
  )

  # Prepare CellChat inputs
  cc_inputs <- get_cc_inputs(
    obj = cond_obj,
    assay = "RNA",
    layer = "data",
    sample_col = SAMPLE_COL,
    group_cols = GROUP_COLS,
    coord_cols = COORD_COLS,
    conversion_factor = CONVERSION_FACTOR
  )

  # Seurat object is no longer needed for this condition
  rm(cond_obj)
  gc()

  # Create CellChat object
  cond_cc_obj <- createCellChat(
    object = cc_inputs$data_input,
    meta = cc_inputs$meta,
    datatype = "spatial",
    group.by = NICHE_COL,
    coordinates = cc_inputs$spatial_locs,
    spatial.factors = cc_inputs$spatial_factors
  )

  # Inputs no longer needed
  rm(cc_inputs)
  gc()

  # Drop unused metadata levels
  cond_cc_obj <- cellchat_droplevels(
    cond_cc_obj,
    meta_cols = GROUP_COLS,
    clear_idents = TRUE
  )

  # Run CellChat
  cond_cc_obj <- cellchat_pipeline(
    cond_cc_obj,
    db = cellchatDB.use,
    contact_range = CONTACT_RANGE,
    interaction_range = INTERACTION_RANGE,
    scale_distance = SCALE_DISTANCE
  )

  # Save
  saveRDS(
    cond_cc_obj,
    file.path(
      OUT_DIR,
      sprintf("%s.cellchat.rds", cond)
    )
  )

  # Free memory before moving to next condition
  rm(cond_cc_obj)
  gc()
}