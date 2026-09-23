
# %% Run Libraries 
library(Seurat)
library(SeuratObject)
library(dplyr)
library(pryr) # to check ram usage mem_used()

# %% setting paths
data_dir <- "/mnt/data/project0062/cosmx_gray/data"

slides_dir <- file.path(data_dir, "slides_rds")
meta_path <- file.path(data_dir, "combined/meta/Combined_v0.csv")
comb_rds_path <- file.path(data_dir, "combined/rds/Combined_v0.RDS")

# %% Load Seurat objects

# list all the files in the directory
slides_files <- list.files(slides_dir, pattern = "*.RDS", full.names = TRUE)

obj_list <- lapply(slides_files, readRDS)

# Merge all Seurat objects

combined_obj <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  project = "Lung"
)

# %% check metadata

glimpse(combined_obj@meta.data)

# get slide run tissue name
meta <- combined_obj@meta.data

# %% assigning fov names

# [LOGIC]
# We do this becuase we have multiple slides
# some fov numbers are the same
meta$fov_names <- paste(meta$slide_ID_numeric, meta$fov, sep="_")

# %% assigning slide names
meta$slide_name <- paste0("slide", meta$slide_ID_numeric)

# %% Rename slide image names

# [INFO]
# Run_Tissue_name column has the same values as the object @image slot names
# but with spaces and slashes replaced by dots

# get all the run tissue names
run_tissue_names <- unique(meta$Run_Tissue_name)

# get all the @images slot names
slide_image_names <- names(combined_obj@images)


# [LOGIC]
# Now we want to rename the @images slot names
# with the slide_names corresponding to the Run_Tissue_name

slide_names_dict <- list()
for (rt_name in run_tissue_names) {
  # get the corresponding slide name
  slide_name <- meta$slide_name[meta$Run_Tissue_name == rt_name][1]
  # remove spaces and slashes from the Run_Tissue_name
  # and replace with dots
  rt_name_clean <- gsub("[/ ]", ".", rt_name)
  # add to the dictionary
  slide_names_dict[[rt_name_clean]] <- slide_name
}

# Get the current names of the @images slot in the Seurat object
old_names <- names(combined_obj@images)

# Loop over each image name
for (i in seq_along(old_names)) {
  # Store the current (old) image name
  old_name <- old_names[i]

  # Check if this image name is in the renaming dictionary
  if (old_name %in% names(slide_names_dict)) {
    # Get the new name from the dictionary
    new_name <- slide_names_dict[[old_name]]
    # Assign the image to the new name in the list
    combined_obj@images[[new_name]] <- combined_obj@images[[old_name]]
    # Update the internal key (important for features to match correctly)
    combined_obj@images[[new_name]]@key <- paste0(new_name, "_")
    # Remove the old image name to complete the rename
    combined_obj@images[[old_name]] <- NULL
  }
}

# %% update the object metadata and export

# update metadata
combined_obj@meta.data <- meta

# export metadata
write.csv(meta, file = meta_path, row.names = TRUE)

# %% Save the merged object as RDS
saveRDS(combined_obj, file = comb_rds_path)