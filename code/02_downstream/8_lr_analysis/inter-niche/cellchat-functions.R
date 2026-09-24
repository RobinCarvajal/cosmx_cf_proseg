
#### Helper functions #####

# Compute sample-specific spatial factors for CellChat
sf_per_sample <- function(meta, spatial_locs, sample_col, conversion_factor) {

  # Standardize coordinate column names so CellChat knows these are x/y positions
  names(spatial_locs) <- c("x", "y")
  
  # Initialize empty data frame to store one row of spatial factors per sample
  spatial_factors <- data.frame()
  
  # Get the unique sample names from the selected metadata column
  sample_names <- unique(meta[[sample_col]])

  # Loop through each sample separately
  for (sample in sample_names) {
    
    # Identify cells/spots that belong to the current sample
    idx <- meta[[sample_col]] == sample
    
    # Subset spatial coordinates to only the current sample
    spatial_locs_sample <- spatial_locs[idx, ]
    
    # Compute pairwise cell-to-cell distances for this sample
    d <- CellChat::computeCellDistance(spatial_locs_sample)
    
    # Use the minimum observed distance as a proxy for nearest-neighbour spacing
    # and convert it from pixels to micrometers using the provided conversion factor
    spot_size <- min(d) * conversion_factor
    
    # Store the spatial parameters for this sample:
    # - ratio: pixel-to-micrometer conversion factor
    # - tol: half the inferred spot size, used as spatial tolerance
    spatial_factors_sample <- data.frame(
      sample = sample, 
      ratio = conversion_factor, 
      tol = spot_size / 2
    )

    print(spatial_factors_sample)

    # Append the current sample result to the full output table
    spatial_factors <- rbind(spatial_factors, spatial_factors_sample)
  }

  # Use sample names as row names for easier downstream matching
  row.names(spatial_factors) <- spatial_factors$sample

  # Return one row of spatial factors per sample
  return(spatial_factors)
}


###### Function that extracts and prepares CellChat inputs ################
get_cc_inputs <- function(
  obj, 
  assay = "RNA", 
  layer = "data", 
  sample_col = "sample_name",
  group_cols = c("condition"),
  coord_cols = c("CenterX_global_px", "CenterY_global_px"),
  conversion_factor = 0.12028
){
  # Extract expression matrix from the selected assay/layer
  # This will be used as the main input matrix for CellChat
  data_input <- GetAssayData(obj, assay = assay, layer = layer)

  # Extract metadata columns needed for CellChat grouping
  # Includes the sample identifier plus any biological grouping variables
  meta <- obj@meta.data[c(sample_col, group_cols)]
  # make all columns factors for downstream consistency
  meta[] <- lapply(meta, factor)

  # Extract spatial coordinates from metadata
  # These coordinates are expected to be in pixel units
  spatial_locs <- obj@meta.data[coord_cols]
  
  # Rename coordinate columns to x/y for downstream consistency
  colnames(spatial_locs) <- c("x", "y")

  # Compute sample-specific spatial factors using the metadata,
  # spatial coordinates, sample labels, and pixel-to-micrometer ratio
  spatial_factors <- sf_per_sample(
    meta, 
    spatial_locs, 
    sample_col, 
    conversion_factor
  )

  # rename the sample column to "sample" for downstream consistency
  colnames(meta)[colnames(meta) == sample_col] <- "samples"

  # Bundle all prepared inputs into a single list for downstream use
  out <- list(
    data_input = data_input,
    meta = meta,
    spatial_locs = spatial_locs,
    spatial_factors = spatial_factors
  )

  # Return all CellChat-ready inputs
  return(out)

}

# cellchat_drop_levels ---------------------------------------------------------

cellchat_droplevels <- function(data, meta_cols, clear_idents = TRUE){
  #parameter# data: CellChat object
  #parameter# meta_cols: columns to drop levels
  #parameter# clear_idents: clear idents
  #output#    cleared_data: cleared CellChat object
  
  # dropping levels of idents
  if (clear_idents == TRUE){
    
    data@idents <- droplevels(data@idents)
    
  }
  
  # dropping levels in a loop 
  for (col in meta_cols){
    data@meta[[col]] <- droplevels(data@meta[[col]])
  }
  
  # return the result as a cleared output
  cleared_data <- data
  
  return(cleared_data)
}

# cellchat_pipeline ------------------------------------------------------------

cellchat_pipeline <- function(
  data, 
  db, 
  contact_range, 
  interaction_range, 
  scale_distance
){
  #parameter# obj: CellChat object
  #parameter# db: CellChat database
  #parameter# contact_range: cell soma/body diameter (mirons)
  #parameter# interaction_range: the maximum interaction/diffusion length of ligands (microns)
  
  # Adding DB
  data@DB <- db
  
  # Subsetting data to include only genes involved in cell interactions
  data <- subsetData(data)
  
  # Standard pipeline
  data <- identifyOverExpressedGenes(data)
  print("Over expressed genes identified")
  data <- identifyOverExpressedInteractions(data)
  print("Over expressed interactions identified")

  
  # compute communication probability  
  data <- computeCommunProb(data, type="truncatedMean",
                            interaction.range=interaction_range, 
                            contact.range=contact_range,
                            scale.distance=scale_distance
                            )
  # the rest of the parameters are set to default 
  
  print("Communication Probability computed")
  data <- filterCommunication(data, min.cells = 10)
  print("Cummunications filtered")
  data <- computeCommunProbPathway(data)
  print("Communication Probability Pathway computed")
  data <- aggregateNet(data)
  print("Aggregated Cell-Cell communication network calculated")
  # latest addition 
  data <- netAnalysis_computeCentrality(data, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  print("Network centrality scores computed")

  
  return(data)
}