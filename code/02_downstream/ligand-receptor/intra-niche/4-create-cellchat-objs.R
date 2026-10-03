# The purpose of this script is to create CellChat objects for each condition 
# in the dataset, which will be used for ligand-receptor interaction analysis.

# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
library(SpatialCellChat)
#library(CellChat)
library(patchwork)
library(future)
options(stringsAsFactors = FALSE)


# Set paths --------------------------------------------------------------------
MAIN_DIR <- "/mnt/data/project0062/proseg_data"
setwd(MAIN_DIR)
source("code/8_lr_analysis/cellchat-functions.R")

# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Set up parallel processing ---------------------------------------------------
# future::plan("multisession", workers = 4) 
# options(future.globals.maxSize = 4 * 1024^3)
# #helps speed up LR identification

# Niche object path ------------------------------------------------------------

IN_OBJ_PATH <- file.path(
  MAIN_DIR, "results", "comb-full", "lr-manual-niches", "airway_epithelium", "niche-obj.rds"
)

OUT_DIR <- file.path(
  MAIN_DIR, "results", "comb-full", "lr-manual-niches", "airway_epitheliu "
)

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
NICHE_COL <- "niche_name"
SAMPLE_COL <- "sample_name"
GROUP_COLS <- c(CONDITION_COL, NICHE_COL, CELLTYPE_COL)
COORD_COLS <- c("centroid_x", "centroid_y")

obj_list <- SplitObject(niche_obj, split.by = CONDITION_COL)


# CellChat inputs --------------------------------------------------------------

## Need
# data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
# meta <- obj@meta.data[group_vars]
# spatial.locs <-  obj@meta.data[coord_vars]
# contact.range <- 10


# Prepare all the inputs -------------------------------------------------------

# There is one object per condition `obj_list`
cond_names <- names(obj_list)

# Initialize an empty list to store processed inputs for each object
obj_input_list <- list()

# Loop through each object by name
for (cond in cond_names) {
  
  # Retrieve the current object from the list
  cond_obj <- obj_list[[cond]]
  
  # Extract and prepare the inputs needed for CellChat
  cc_inputs <- get_cc_inputs(
    obj=cond_obj, 
    assay = "RNA", 
    layer = "data", 
    sample_col = SAMPLE_COL,
    group_cols = GROUP_COLS,
    coord_cols = COORD_COLS,
    conversion_factor = CONVERSION_FACTOR
  )
  
  # Store the prepared inputs in the output list using the same object name
  obj_input_list[[cond]] <- cc_inputs
}

# Create the cellchat objects --------------------------------------------------

# This will create one CellChat object per condition
# using the prepared inputs from the previous step

# Loop through each set of inputs by object name
for (cond in cond_names) {
  # Retrieve the prepared inputs for the current object
  cc_inputs <- obj_input_list[[cond]]
  
  # Create a CellChat object using the prepared inputs
  cond_cc_obj <- createSpatialCellChat(
    object = cc_inputs$data_input,
    meta = cc_inputs$meta,
    datatype = "spatial",
    group.by = CELLTYPE_COL,
    coordinates = cc_inputs$spatial_locs,
    spatial.factors = cc_inputs$spatial_factors,
  )
  message()
  
  # dropping unnecessary levels from the metadata
  clean_cond_cc_obj <- cellchat_droplevels(
    cond_cc_obj, 
    meta_cols = GROUP_COLS,
    clear_idents = TRUE
  )

  # run cellchat pipeline for interactions
  processed_cond_cc_obj <- cellchat_pipeline(
    clean_cond_cc_obj, 
    db = cellchatDB.use,
    contact_range = CONTACT_RANGE, 
    interaction_range = INTERACTION_RANGE,
    scale_distance = SCALE_DISTANCE
  )

  break
  # Save the cond cc object
  saveRDS(processed_cond_cc_obj, file.path(OUT_DIR, sprintf("%s.cellchat.rds", cond)))
}

