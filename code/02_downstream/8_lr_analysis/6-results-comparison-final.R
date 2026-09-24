# Set Working Dir --------------------------------------------------------------
MAIN_DIR <- "/data/cosmx_gray"
setwd(MAIN_DIR)

# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
#library(SpatialCellChat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
source("code/8_lr_analysis/cellchat-functions.R")

# call reticulate from conda environment with umap learn isntalled
library(reticulate)

use_condaenv(
  "cellchat-env",
  conda = "/data/miniforge3/bin/conda",
  required = TRUE
)

py_config()

# Set paths --------------------------------------------------------------------

# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Niche object path ------------------------------------------------------------
IN_DIR <- file.path(
  MAIN_DIR, "results", "comb", "lr-manual-niches", "1"
)

OUT_DIR <- file.path(
  MAIN_DIR, "results", "comb", "lr-manual-niches", "1", "CF_vs_CF_ETI"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Define input CellChat objects ------------------------------------------------
# Adjust these names if needed
obj_paths <- c(
  CTRL = file.path(IN_DIR, "CF.cellchat.rds"),
  CF   = file.path(IN_DIR, "CF_ETI.cellchat.rds")
)

# Check files exist
missing_paths <- obj_paths[!file.exists(obj_paths)]
if (length(missing_paths) > 0) {
  stop(
    "The following CellChat object(s) were not found:\n",
    paste(names(missing_paths), missing_paths, sep = ": ", collapse = "\n")
  )
}

# Load CellChat objects --------------------------------------------------------
message("Loading CellChat objects...")
cellchat_list <- lapply(obj_paths, readRDS)

# Ensure names are preserved
names(cellchat_list) <- names(obj_paths)

# get all the labels
group_new <- union(
  levels(cellchat_list[[1]]@idents), levels(cellchat_list[[2]]@idents)
)

# lift the cellchat objects to have the same labels
cellchat_list <- lapply(cellchat_list, function(x) {
  liftCellChat(x, group.new = group_new)
})

### Sanity Check ###

celltypes_1 <- rownames(cellchat_list[[1]]@net$count)
celltypes_2 <- rownames(cellchat_list[[2]]@net$count)

if (!identical(celltypes_1, celltypes_2)) {
  stop(
    paste0(
      "Skipping netVisual_diffInteraction: cell types do not match.\n",
      "Only in object 1: ", paste(setdiff(celltypes_1, celltypes_2), collapse = ", "), "\n",
      "Only in object 2: ", paste(setdiff(celltypes_2, celltypes_1), collapse = ", ")
    )
  )
}


#######
# Compute centrality scores for each CellChat object (if not already done) -----
for (item in names(cellchat_list)) {
  # calculate centrality scores if not already done
  cellchat_list[[item]] <- netAnalysis_computeCentrality(
    cellchat_list[[item]],
    slot.name = "netP"
  )
}


# Merge CellChat objects -------------------------------------------------------
message("Merging CellChat objects...")
cellchat_merged <- mergeCellChat(
  object.list = cellchat_list,
  add.names = names(cellchat_list)
)

#### PART 1: IDENTIFY ALTERED INTERACTIONS AND CELL POPULATIONS ####

# 1.1. compare the total number of interactions and interaction stregth

gg1 <- compareInteractions(
  cellchat_merged, 
  show.legend = FALSE, 
  group = c(1, 2),
  measure = "count"
)
gg2 <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

png(
  filename = file.path(OUT_DIR, "barplot-diff-summary.png"),
  width = 1800, height = 900, res = 300
)
gg1+gg2
dev.off()

# 1.2. compare the number of interactions and interaction strength among different cell populations

# Red = increased in second dataset in group comparison
# Blue = decreased in second dataset in group comparison

# Compare Interactions Networks - Count
png(
  filename = file.path(OUT_DIR, "circle-diff-count.png"),
  width = 1600, height = 1600, res = 300
)

par(xpd=TRUE)

netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "count",
  label.edge = FALSE
)

dev.off()

# Compare Interactions Networks - Weight

# open device
png(
  filename = file.path(OUT_DIR, "circle-diff-strength.png"),
  width = 1600, height = 1600, res = 300
)
# plot
par(xpd=TRUE)
netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight",
  label.edge = FALSE
)
# close device
dev.off()

# Heatmaps of differential interactions ---------------------------------------

# joint plot for strength and interaction

# open device
png(
  filename = file.path(OUT_DIR, "heatmap-diff-count-strength.png"),
  width =15, height = 7, res = 300, units = "in"
)
# count plot
gg1 <- netVisual_heatmap(cellchat_merged, measure = "count")
# strength plot
gg2 <- netVisual_heatmap(cellchat_merged, measure = "weight")
# combine and plot
gg1 + gg2
# close device
dev.off()


# Compare the number of interactions and interaction strength among different 
# cell populations

# [C] Circle plot showing the number of interactions or interaction
#  strength among different cell populations across multiple datasets


# Count

# open device
png(
  file.path(OUT_DIR, "circle-comp-count.png"),
  width = 2400, height = 1200, res = 180
)
# plot
weight.max <- getMaxWeight(cellchat_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(
    cellchat_list[[i]]@net$count, 
    weight.scale = T, 
    label.edge= F, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(cellchat_list)[i])
  )
}
# close device
dev.off()

# Stength
# open device
png(
  file.path(OUT_DIR, "circle-comp-strength.png"),
  width = 2400, height = 1200, res = 180
)
# plot
weight.max <- getMaxWeight(cellchat_list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(
    cellchat_list[[i]]@net$weight, 
    weight.scale = T, 
    label.edge= F, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Interaction strength - ", names(cellchat_list)[i])
  )
}
# close device
dev.off()

# D) Circle plot showing the differntial number of interactions or interaction 
# strength among different cell populations across multiple datasets

# define the coarse cell type groups
group.cellType <- c(
  "MYELOID",      # alveolar_macrophages
  "OTHER",        # ambiguous
  "EPI",   # AT1
  "EPI",   # AT2
  "LYMPHOID",     # B
  "EPI",   # basal
  "EPI",   # ciliated
  "ENDOTHELIAL",  # endothelial
  "STROMAL",      # fibroblasts
  "MYELOID",      # macrophages
  "MYELOID",      # mast
  "MYELOID",      # neutrophils
  "LYMPHOID",     # plasma
  "EPI",   # secretory
  "EPI",   # secretory-ciliated
  "STROMAL",      # smooth_muscle
  "OTHER",        # sputum
  "LYMPHOID"      # T
)

# set the factor levels to control the order in the plot
group.cellType <- factor(
  group.cellType,
  levels = c("EPI", "MYELOID", "LYMPHOID", "ENDOTHELIAL", "STROMAL", "OTHER")
)

# merge the interactions based on the defined cell type groups
cellchat_list_grouped <- lapply(cellchat_list, function(x) {
  mergeInteractions(x, group.cellType)
})

# merge the grouped CellChat objects for later use
cellchat_grouped <- mergeCellChat(
  cellchat_list_grouped,
  add.names = names(cellchat_list_grouped)
)

# control the edge weights across different datasets
weight.max <- getMaxWeight(
  cellchat_list_grouped,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "count", "count.merged")
)

# open device
png(
  file.path(OUT_DIR, "circle-grouped-comp-count.png"),
  width = 2400, height = 1200, res = 180
)
# plot
par(mfrow = c(1, length(cellchat_list_grouped)), xpd = TRUE)
for (i in seq_along(cellchat_list_grouped)) {
  netVisual_circle(
    cellchat_list_grouped[[i]]@net$count.merged,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = weight.max[3],
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", names(cellchat_list_grouped)[i]),
    margin = 0.5
  )
}
# close device
dev.off()

# control the edge weights across different datasets
weight.max.str <- getMaxWeight(
  cellchat_list_grouped,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "weight", "weight.merged")
)

# open device
png(
  file.path(OUT_DIR, "circle-grouped-comp-strength.png"),
  width = 2400, height = 1200, res = 180
)
# plot
par(mfrow = c(1, length(cellchat_list_grouped)), xpd = TRUE)

for (i in seq_along(cellchat_list_grouped)) {
  netVisual_circle(
    cellchat_list_grouped[[i]]@net$weight.merged,
    weight.scale = TRUE,
    label.edge = FALSE, # not very informative as values are relatively small after merging
    edge.weight.max = weight.max.str[3],
    edge.width.max = 12,
    title.name = paste0("Interaction strength - ", names(cellchat_list_grouped)[i]),
    margin = 0.5
  )
}

dev.off()


# [GROUPED] Differnetial number of interactions and interaction strength among different 
# cell populations across multiple datasets

# open device
png(
  file.path(OUT_DIR, "circle-grouped-diff-count.png"),
  width = 10, height = 10, res = 300, units = "in"
) 
# count plot
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(
  cellchat_grouped,
  weight.scale = TRUE,
  measure = "count.merged",
  label.edge = FALSE
)
# close device
dev.off()

# open device
png(
  file.path(OUT_DIR, "circle-grouped-diff-strength.png"),
  width = 10, height = 10, res = 300, units = "in"
)
# strength plot
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(
  cellchat_grouped,
  weight.scale = TRUE,
  measure = "weight.merged",
  label.edge = FALSE
)
# close device
dev.off()


## 1.3) COMPARE THE MAJOR SOURCES AND TARGETS IN A 2D SPACE ##

### [A] Identify cell populations with significant changes in sending or
### receiving signals 

# total interactions per cell group (in + out, minus self)
num.link <- sapply(
  cellchat_list, 
  function(x) {
    rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)
  }
)
# global min/max to standardize dot sizes across datasets
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
# list to store plots
gg <- list()
# generate signaling role scatter plots per dataset
for (i in 1:length(cellchat_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    cellchat_list[[i]], 
    title = names(cellchat_list)[i], 
    weight.MinMax = weight.MinMax
  )
}
# combine and plot
plot <- patchwork::wrap_plots(plots = gg)

# open device
png(
  filename = file.path(OUT_DIR, "scatter-comp-signalling-role.png"),
  width = 10, height = 4, res = 300, units = "in"
)
# plot
print(plot)
# close device
dev.off()

### [B] Identify the signaling changes of specific cell populations

# get the cell types from the merged object
celltypes <- levels(cellchat_merged@meta$ann_lvl_3)

# DIR for the scatter plots of signaling changes for each cell type
signalling_changes_dir <- file.path(OUT_DIR, "scatter-signalling-changes")
dir.create(signalling_changes_dir, showWarnings = FALSE)

# loop through cell types and generate scatter plots for each, with error handling
for (celltype in celltypes) {

  # if fails next celltype
  tryCatch({
    gg <- netAnalysis_signalingChanges_scatter(
      cellchat_merged, idents.use = celltype
    )

    png(
      filename = file.path(
        signalling_changes_dir, paste0("scatter-signalling-changes-", celltype, ".png")
      ),
      width = 10, height = 8, res = 180, units = "in"
    )
    print(gg)
    # close device
    dev.off()

  }, error = function(e) {
    message("Failed for celltype: ", celltype)

    # save the failed celltype in a log file
    log_file <- file.path(signalling_changes_dir, "scatter_signalling_changes_failures.log")
    write(
      paste("Failed for celltype:", celltype), 
      file = log_file, 
      append = TRUE
    )
    
  })
}


################################################################################
# PART 2: IDENTIFY ALTERED SINGALLING WITH DISTINCT NETWORK ARCHITECTURE AND   # 
# INTERACTION STRENGTH                                                         #
################################################################################

## 2.1) IDENTIFY SIGNALLING NETWORKS WITH LARGER (OR LESS) DIFFERENCE AS WELL AS 
## SIGNALLING GROUPS BASED ON THEIR FUNCTIONAL/STRUCTURE SIMILARITY

### Identify signalling groups based on their functional similarity

# Compute signaling network similarity for any pair of datasets
cellchat_merged <- computeNetSimilarityPairwise(
  cellchat_merged, type = "functional"
)
# Manifold learning of the signaling networks based on their similarity
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
# Classification learning of the signaling networks
cellchat_merged <- netClustering(cellchat_merged, type = "functional")

# DIR for storing the results of manifold learning and similarity analysis
umap_signalling_similarity_dir <- file.path(
  OUT_DIR, "umap_signalling_similarity"
)
dir.create(umap_signalling_similarity_dir, showWarnings = FALSE)


# open device
png(
  filename = file.path(
    umap_signalling_similarity_dir, "umap-similarity-functional.png"
  ),
  width = 1800, height = 1400, res = 180
)
# 2D visualization of the joint manifold learning of signaling networks 
# from two datasets
netVisual_embeddingPairwise(
  cellchat_merged, type = "functional", label.size = 3.5
)
# close device
dev.off()

# open device
png(
  filename = file.path(
    umap_signalling_similarity_dir, "umap-zoom-similarity-functional.png"
  ),
  width = 1800, height = 1400, res = 180
)
# Zoom-in visualization of the joint manifold learning of signaling networks
netVisual_embeddingPairwiseZoomIn(
  cellchat_merged, type = "functional", nCol = 2
)
# close device
dev.off()

### Identify signalling groups based on structure similarity

# Compute signaling network similarity for any pair of datasets
cellchat_merged <- computeNetSimilarityPairwise(
  cellchat_merged, type = "structural"
)
# Manifold learning of the signaling networks based on their similarity
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
# Classification learning of the signaling networks
cellchat_merged <- netClustering(cellchat_merged, type = "structural")

# open device
png(
  filename = file.path(
    umap_signalling_similarity_dir, "umap-similarity-structural.png"
  ),
  width = 1800, height = 1400, res = 180
)
# Visualization in 2D-space
netVisual_embeddingPairwise(
  cellchat_merged, type = "structural", label.size = 3.5
)
# close device
dev.off()

png(
  filename = file.path(
    umap_signalling_similarity_dir, "umap-zoom-similarity-structural.png"
  ),
  width = 1800, height = 1400, res = 180
)
# Zoom-in visualization of the joint manifold learning of signaling networks
netVisual_embeddingPairwiseZoomIn(
  cellchat_merged, type = "structural", nCol = 2
)
# close device
dev.off()


### Compute and visualize the pathway distance in the learned joint manifold

# Rank the similarity of the shared signaling pathways based on their joint 
# manifold learning

# functional 

# open device
png(
  filename = file.path(
    umap_signalling_similarity_dir, "rank-similarity-functional.png"
  ),
  width = 1800, height = 1400, res = 180
)
# plot 
rankSimilarity(
  cellchat_merged, 
  type = "functional",
  title = "Functional similarity of signaling pathways"
)
# close device
dev.off()

# structural

# open device
png(
  filename = file.path(
    umap_signalling_similarity_dir, "rank-similarity-structural.png"
  ),
  width = 1800, height = 1400, res = 180
)
# plot
rankSimilarity(
  cellchat_merged, 
  type = "structural",
  title = "Structural similarity of signaling pathways"
)
# close device
dev.off()

## 2.2) Identify altered signalling with distinct interaction strength

### [A] Compare the overall information flow of each signaling pathway or 
### ligand receptor pair

gg1 <- rankNet(
  cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, 
  targets.use = NULL, stacked = T, do.stat = TRUE
)
gg2 <- rankNet(
  cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, 
  targets.use = NULL, stacked = F, do.stat = TRUE
)

# open device
png(
  filename = file.path(OUT_DIR, "barplot-pathway-contribution.png"),
  width = 10, height = 5, res = 300, units = "in"
)
# plot
gg1 + gg2
# close device
dev.off()

### [B] Compare outgoing (or incoming) signalling pattern associated with
### each cell population


library(ComplexHeatmap)

# numeric index for the datasets to compare
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(cellchat_list[[i]]@netP$pathways, cellchat_list[[i+1]]@netP$pathways)

# small fix so that the heatmap plotting function does not crash
#
# explanation:
# when using the union of pathways across datasets, some pathways may be
# absent in a given CellChat object, resulting in matrices with identical
# values (commonly all zeros). Internally, the heatmap function generates
# a color scale using the minimum and maximum values, but this fails when
# both are identical:
#
#   circlize::colorRamp2(c(0, 0), ...)
#
# adding a negligible pseudocount introduces minimal variation without
# meaningfully affecting the biological interpretation or visualization

epsilon <- 1e-10

for (pathway in names(cellchat_list[[i]]@netP$centr)) {

  cellchat_list[[i]]@netP$centr[[pathway]]$outdeg <-
    cellchat_list[[i]]@netP$centr[[pathway]]$outdeg + epsilon
}

# outgoing signaling pattern

# open device
png(
  filename = file.path(OUT_DIR, "heatmap-signalling-patterns-outgoing.png"),
  width = 7, height = 7, res = 300, units = "in"
)
# plots
ht1 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i]], pattern = "outgoing", signaling = pathway.union, 
  title = names(cellchat_list)[i], width = 5, height = 6
)
ht2 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i+1]], pattern = "outgoing", signaling = pathway.union, 
  title = names(cellchat_list)[i+1], width = 5, height = 6
)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# close device
dev.off()

# incoming signaling pattern

# open device
png(
  filename = file.path(OUT_DIR, "heatmap-signalling-patterns-incoming.png"),
  width = 15, height = 10, res = 180, units = "in"
)
# plots
ht1 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i]], pattern = "incoming", signaling = pathway.union, 
  title = names(cellchat_list)[i], width = 5, height = 6
)
ht2 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i+1]], pattern = "incoming", signaling = pathway.union, 
  title = names(cellchat_list)[i+1], width = 5, height = 6
)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# close device
dev.off()

# all signalling patterns

# open device
png(
  filename = file.path(OUT_DIR, "heatmap-signalling-patterns-all.png"),
  width = 15, height = 10, res = 180, units = "in"
)
# plots
ht1 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i]], pattern = "all", signaling = pathway.union, 
  title = names(cellchat_list)[i], width = 5, height = 6, color.heatmap = "OrRd"
)
ht2 = netAnalysis_signalingRole_heatmap(
  cellchat_list[[i+1]], pattern = "all", signaling = pathway.union, 
  title = names(cellchat_list)[i+1], width = 5, height = 6, color.heatmap = "OrRd"
)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# close device
dev.off()


###################################################################
# PART 3: IDENTIFY THE UP-REGULATED AND DOWN-REGULATED SIGNALLING #
# LIGAND-RECEPTOR PAIRS                                           # 
###################################################################


## 3.1) Identify dysfunctional signaling by comparing the communication
## probabilities


# DIR for storing the bubble plots of up and down regulated interactions
bubble_signaling_changes_strength_dir <- file.path(
  OUT_DIR, "bubble-signaling-changes-strength"
)
dir.create(bubble_signaling_changes_strength_dir, showWarnings = FALSE)

# loop through cell types and generate bubble plots for each
for (i in 1:length(celltypes)) {

    tryCatch({

      # increased signalling in CF
      gg <- netVisual_bubble(
        cellchat_merged, 
        sources.use = i, 
        comparison = c(1, 2), 
        max.dataset = 2, 
        title.name = paste0("Increased signaling in CF - ", celltypes[i]), 
        angle.x = 45, 
        remove.isolate = TRUE
      )
    
      # open device
      png(
        filename = file.path(
          bubble_signaling_changes_strength_dir, 
          paste0("bubble-signaling-changes-up-", celltypes[i], ".png")
        ),
        width = 10, height = 8, res = 180, units = "in"
      )
      # plot
      print(gg)
      # close device
      dev.off()

    }, error = function(e) {
      message("Failed for celltype: ", celltypes[i])

      # save the failed celltype in a log file
      log_file <- file.path(
        bubble_signaling_changes_strength_dir, 
        "failed-up.log"
      )
      write(
        paste("Failed for celltype:", celltypes[i]), 
        file = log_file, 
        append = TRUE
      )

    })
}


# loop through cell types and generate bubble plots for each
for (i in 1:length(celltypes)) {

    tryCatch({
    
      # decreased signalling in CF
      gg <- netVisual_bubble(
        cellchat_merged, 
        sources.use = i, 
        comparison = c(1, 2), 
        max.dataset = 1, 
        title.name = paste0("Decreased signaling in CF - ", celltypes[i]), 
        angle.x = 45, 
        remove.isolate = TRUE
      )
    
      # open device
      png(
        filename = file.path(
          bubble_signaling_changes_strength_dir, 
          paste0("bubble-signaling-changes-down-", celltypes[i], ".png")
        ),
        width = 10, height = 8, res = 180, units = "in"
      )
      # plot
      print(gg)
      # close device
      dev.off()

    }, error = function(e) {
      message("Failed for celltype: ", celltypes[i])

      # save the failed celltype in a log file
      log_file <- file.path(
        bubble_signaling_changes_strength_dir, 
        "failed-down.log"
      )
      write(
        paste("Failed for celltype:", celltypes[i]), 
        file = log_file, 
        append = TRUE
      )

    })
}


### 3.2) Identify dysfunctional signaling by using diff. exp. analysis


# we load the cellchat_list again and subset to include cell types 
# with enough cells in both datasets for performing diff. exp. analysis,
# as the previous merged object is subsetted for visualizing the signaling 
# changes with distinct interaction strength, which may not be suitable 
# for diff. exp. analysis.

message("Loading CellChat objects...")
cellchat_list_de <- lapply(obj_paths, readRDS)

# Ensure names are preserved
names(cellchat_list_de) <- names(obj_paths)


# get labels to preserve
tab_joint <- table(
  cellchat_merged@idents$joint,
  cellchat_merged@meta$datasets
)
# keep cell types with at least 10 cells in both datasets
keep_joint <- rownames(tab_joint)[
  tab_joint[, "CTRL"] >= 10 &
  tab_joint[, "CF"] >= 10
]


# define a positive dataset, i.e., the dataset with positive fold change 
# against the other dataset
pos_dataset = "CF"
# define a char name used for storing the results of diff. exp. analysis
features_name = paste0(pos_dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an 
# ultra-fast Wilcoxon test using the presto package, which gives smaller 
# values of logFC. Thus we here set a smaller value of thresh.fc compared 
# to the original one (thresh.fc = 0.1). Users can also provide a vector and 
# dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged 
# and cellchat@var.features$LS.merged.info. 



for (i in 1:length(cellchat_list_de)) {
  cellchat_list_de[[i]] <- subsetCellChat(
    cellchat_list_de[[i]],
    idents.use = keep_joint
  )
}

cellchat_merged_de <- mergeCellChat(
  object.list = cellchat_list_de,
  add.names = names(cellchat_list_de)
)


cellchat_merged_de <- identifyOverExpressedGenes(
  cellchat_merged_de, 
  group.dataset = "datasets", 
  pos.dataset = pos_dataset, 
  features.name = features_name, 
  only.pos = FALSE, 
  thresh.pc = 0.1, 
  thresh.fc = 0.05,
  thresh.p = 0.05, 
  group.DE.combined = FALSE
)

# map the results of differential expression analysis onto the inferred 
# cell-cell communications to easily manage/subset the ligand-receptor 
# pairs of interest
net <- netMappingDEG(
  cellchat_merged_de, 
  features.name = features_name,
  variable.all = TRUE
)

# extract the ligand-receptor pairs with upregulated ligands in CF
net_up <- subsetCommunication(
  cellchat_merged_de, 
  net = net, 
  datasets = "CF",
  ligand.logFC = 0.05, 
  receptor.logFC = NULL
)

# extract the ligand-receptor pairs with upregulated ligands and upregulated 
# receptors in CTRL, i.e.,downregulated in CF
net_down <- subsetCommunication(
  cellchat_merged_de, 
  net = net, 
  datasets = "CTRL",
  ligand.logFC = -0.05, 
  receptor.logFC = NULL
)

# gene lists
gene_up <- extractGeneSubsetFromPair(net_up, cellchat_merged_de)
gene_down <- extractGeneSubsetFromPair(net_down, cellchat_merged_de)


### 3.3) Visualize the identified up-regulated and down-regulated signaling 
### ligand-receptor pairs

# [A] Bubble plot

pairLR_use_up = net_up[, "interaction_name", drop = F]
sources_up <- unique(net_up$source)
sources_up_idx <- which(celltypes %in% sources_up)


# DIR for diff. exp. bubble plots
bubble_signaling_changes_de_dir <- file.path(
  OUT_DIR, "bubble-signaling-changes-de"
)
dir.create(bubble_signaling_changes_de_dir, showWarnings = FALSE)


# upregulated signalling

for (i in sources_up_idx) {

  tryCatch({
    
    message("Processing source cell type: ", celltypes[i])

    # subset the net_up for the current source cell type
    net_up_celltype <- net_up[net_up$source == celltypes[i], ]
    pairLR_use_up <- net_up_celltype[, "interaction_name", drop = F]
    targets_up <- unique(net_up_celltype$target)
    targets_up_idx <- which(celltypes %in% targets_up)

    gg <- netVisual_bubble(
      cellchat_merged_de, 
      pairLR.use = pairLR_use_up, 
      sources.use = i, 
      targets.use = targets_up_idx,
      comparison = c(1, 2),  
      angle.x = 90, 
      remove.isolate = TRUE,
      title.name = paste0(
        "Up-regulated signaling in ", 
        names(cellchat_list_de)[2], " - Source: ", celltypes[i]
      )
    )
    # open device
    png(
      filename = file.path(
        bubble_signaling_changes_de_dir, 
        paste0("bubble-signaling-changes-up-", celltypes[i], ".png")
      ),
      width = 10, height = 8, res = 180, units = "in"
    )
    # plot
    print(gg)
    # close device
    dev.off()
  }, error = function(e) {
    message("Failed for source cell type: ", celltypes[i])
    
    # save the failed cell type in a log file
    log_file <- file.path(
      bubble_signaling_changes_de_dir, 
      "failed-up.log"
    )
    write(
      paste("Failed for source cell type:", celltypes[i]), 
      file = log_file, 
      append = TRUE
    )
  }
  )
}

# downregulated signalling

pairLR_use_down = net_down[, "interaction_name", drop = F]
sources_down <- unique(net_down$source)
sources_down_idx <- which(celltypes %in% sources_down)
targets_down <- unique(net_down$target)
targets_down_idx <- which(celltypes %in% targets_down)


for (i in sources_down_idx) {

  tryCatch({

    message("Processing source cell type: ", celltypes[i])

    # subset the net_down for the current source cell type
    net_down_celltype <- net_down[net_down$source == celltypes[i], ]
    pairLR_use_down <- net_down_celltype[, "interaction_name", drop = FALSE]
    targets_down <- unique(net_down_celltype$target)
    targets_down_idx <- which(celltypes %in% targets_down)

    gg <- netVisual_bubble(
      cellchat_merged,
      pairLR.use = pairLR_use_down,
      sources.use = i,
      targets.use = targets_down_idx,
      comparison = c(1, 2),
      angle.x = 90,
      remove.isolate = TRUE,
      title.name = paste0(
        "Down-regulated signaling in ",
        names(cellchat_list_de)[2], " - Source: ", celltypes[i]
      )
    )

    png(
      filename = file.path(
        bubble_signaling_changes_de_dir,
        paste0("bubble-signaling-changes-down-", celltypes[i], ".png")
      ),
      width = 10, height = 8, res = 180, units = "in"
    )

    print(gg)
    dev.off()

  }, error = function(e) {

    message("Failed for source cell type: ", celltypes[i])

    log_file <- file.path(
      bubble_signaling_changes_de_dir,
      "failed-down.log"
    )

    write(
      paste("Failed for source cell type:", celltypes[i]),
      file = log_file,
      append = TRUE
    )
  })
}

# NOTE:
# net_up and net_down have the infromation of sender and receptor and lr pair
# it just need to be included in the function, 
# which turns out is a general function for subsetting the 
# communication network based on any criteria, which is very useful 
# for downstream analysis.

# [B] Chord diagram - looks horrible, doesnt give much infromation 

# DIR for chord diagrams of up and down regulated interactions
chord_signalling_changes_de_dir <- file.path(
  OUT_DIR, "chord-signaling-changes-de"
)
dir.create(chord_signalling_changes_de_dir, showWarnings = FALSE)

# upregulated signalling
for (i in sources_up_idx) {
  message("Processing source cell type: ", celltypes[i])
  
  tryCatch({

    # open device
    png(
      filename = file.path(
        chord_signalling_changes_de_dir, 
        paste0("chord-signaling-changes-up-", celltypes[i], ".png")
      ),
      width = 10, height = 10, res = 180, units = "in"
    )
    # plot 
    par(xpd=TRUE)
    netVisual_chord_gene(
      cellchat_list_de[[2]],
      sources.use = i,
      slot.name = 'net',
      net = net_up,
      lab.cex = 0.8,
      small.gap = 3.5,
      title.name = paste0(
        "Up-regulated signaling in ", 
        names(cellchat_list_de)[2], " - Source: ", celltypes[i]
      )
    )
    # close device
    dev.off()

  }, error = function(e) {
    message("Failed for source cell type: ", celltypes[i])
    
    # save the failed cell type in a log file
    log_file <- file.path(
      chord_signalling_changes_de_dir, 
      "failed-up.log"
    )
    write(
      paste("Failed for source cell type:", celltypes[i]), 
      file = log_file, 
      append = TRUE
    )
  }
  )

}

# downregulated signalling

for (i in sources_down_idx) {
  message("Processing source cell type: ", celltypes[i])

  tryCatch({

    # open device
    png(
      filename = file.path(
        chord_signalling_changes_de_dir, 
        paste0("chord-signaling-changes-down-", celltypes[i], ".png")
      ),
      width = 10, height = 10, res = 180, units = "in"
    )
    # plot
    par(xpd=TRUE)
    netVisual_chord_gene(
      cellchat_list_de[[1]],
      sources.use = i,
      slot.name = 'net',
      net = net_down,
      lab.cex = 0.8,
      small.gap = 3.5,
      title.name = paste0(
        "Down-regulated signaling in ", 
        names(cellchat_list_de)[2], " - Source: ", celltypes[i]
      )
    )
    # close device
    dev.off()

  }, error = function(e) {

    message("Failed for source cell type: ", celltypes[i])

    log_file <- file.path(
      chord_signalling_changes_de_dir,
      "failed-down.log"
    )

    write(
      paste("Failed for source cell type:", celltypes[i]),
      file = log_file,
      append = TRUE
    )
  })
}


# [C] Worldcloud plot

# this wordcloud is not very informative, 
# but it can be used to quickly check the most frequent ligands/receptors 
# in the up and down regulated interactions.

# DIR for wordclouds of up and down regulated interactions
wordcloud_signalling_changes_de_dir <- file.path(
  OUT_DIR, "wordcloud-signaling-changes-de"
)
dir.create(wordcloud_signalling_changes_de_dir, showWarnings = FALSE)

# upregulated signalling

# open device
png(
  filename = file.path(
    wordcloud_signalling_changes_de_dir, 
    "wordcloud-signaling-changes-up.png"
  ),
  width = 10, height = 8, res = 180, units = "in"
)
# plot
computeEnrichmentScore(net_up, species = 'human', variable.both = TRUE)
# close device
dev.off()


# downregulated signalling

# open device
png(
  filename = file.path(
    wordcloud_signalling_changes_de_dir, 
    "wordcloud-signaling-changes-down.png"
  ),
  width = 10, height = 8, res = 180, units = "in"
)
# plot
computeEnrichmentScore(net_down, species = 'human', variable.both = TRUE)
# close device
dev.off()


#########################################################################
# PART 4: VISUALLY COMPARE CELL-CELL COMMUNICATION USING HIERACHY PLOT, #
# CIRCLE PLOT, OR CHORD DIAGRAM                                         #
#########################################################################

# get shared pathways across datasets for visualization
pathways_each_dataset <- lapply(cellchat_list, function(x) {
  x@netP$pathways
})
shared_pathways <- Reduce(intersect, pathways_each_dataset)

# visually compare cell-cell communication

# NOTE:
# A pathway needs to be present in both datasets (i.e. conditions) to work

## 4.1) Circle plot 

# DIR for the circle plots for dysregulated pathways
circle_comp_signalling_changes_de_dir <- file.path(
  OUT_DIR, "circle-comp-signaling-changes-de"
)
dir.create(circle_comp_signalling_changes_de_dir, showWarnings = FALSE)

# LOG file for failed visualizations
log_file <- file.path(
  circle_comp_signalling_changes_de_dir, 
  "failed-circle-plot.log"
)
# write header to the log file
write(
  "NOTE: A pathway needs to be present in both datasets (i.e. conditions) to work", 
  file = log_file
)

for (pathway in shared_pathways) {

  pathways.show <- c(pathway)
  print(paste0("Visualizing pathway: ", pathway))

  # loop 
  for (i in 1:length(cellchat_list)) {

    tryCatch({
      
      weight.max <- getMaxWeight(
        cellchat_list, 
        slot.name = c("netP"), 
        attribute = pathways.show
      ) # control the edge weights across different datasets

      # open device 
      png(
        filename = file.path(
          circle_comp_signalling_changes_de_dir, 
          paste0("circle-comp-",pathway,".png")
        ),
        width = 20, height = 10, res = 300, units = "in"
      )

      # plot
      par(mfrow = c(1,2), xpd=TRUE)
      for (i in 1:length(cellchat_list)) {

        netVisual_aggregate(
          cellchat_list[[i]], 
          signaling = pathways.show, 
          layout = "circle", 
          edge.weight.max = weight.max[1], 
          edge.width.max = 10, 
          signaling.name = paste(pathways.show, names(cellchat_list)[i])
        )
      }

      # close device
      dev.off()
      
    }, error = function(e) {
      message(
        "Failed to visualize pathway: ", 
        pathway, 
        " in dataset: ", 
        names(cellchat_list)[i]
      )

      # write to the log file 
      write(
        paste(
          "Failed to visualize pathway:", 
          pathway, 
          "in dataset:", 
          names(cellchat_list)[i]
        ),
        file = log_file,
        append = TRUE
      )
    })

  }
}


## 4.2) Heatmap plot

# DIR for the heatmaps for dysregulated pathways
heatmap_comp_signalling_changes_de_dir <- file.path(
  OUT_DIR, "heatmap-comp-signaling-changes-de"
)
dir.create(heatmap_comp_signalling_changes_de_dir, showWarnings = FALSE)

# LOG file for failed visualizations
log_file <- file.path(
  heatmap_comp_signalling_changes_de_dir, 
  "failed-heatmap-plot.log"
)
# write header to the log file
write(
  "NOTE: A pathway needs to be present in both datasets (i.e. conditions) to work", 
  file = log_file
)

for (pathway in shared_pathways) {

  pathways.show <- c(pathway)
  print(paste0("Visualizing pathway: ", pathway))

  # loop 
  for (i in 1:length(cellchat_list)) {

    tryCatch({
      
      # open device 
      png(
        filename = file.path(
          heatmap_comp_signalling_changes_de_dir, 
          paste0("heatmap-comp-",pathway,".png")
        ),
        width = 20, height = 10, res = 300, units = "in"
      )

      # plot
      par(mfrow = c(1,2), xpd=TRUE)
      ht <- list()
      for (i in 1:length(cellchat_list)) {
        ht[[i]] <- netVisual_heatmap(
          cellchat_list[[i]], 
          signaling = pathways.show, 
          color.heatmap = "Reds",
          title.name = paste(
            pathways.show, "signaling ",names(cellchat_list)[i]
          )
        )
      }
      ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

      # close device
      dev.off()
      
    }, error = function(e) {
      
      msg <- paste(
        "Failed to visualize pathway:", 
        pathway, 
        "in dataset:", 
        names(cellchat_list)[i]
      )
      
      message(msg)
      

      # write to the log file 
      write(
        msg,
        file = log_file,
        append = TRUE
      )
    })

  }
}


## 4.3) Chord diagram

# DIR for the chord diagrams for dysregulated pathways
chord_comp_signalling_changes_de_dir <- file.path(
  OUT_DIR, "chord-comp-signaling-changes-de"
)
dir.create(chord_comp_signalling_changes_de_dir, showWarnings = FALSE)

# LOG file for failed visualizations
log_file <- file.path(
  chord_comp_signalling_changes_de_dir, 
  "failed-chord-plot.log"
)
# write header to the log file
write(
  "NOTE: A pathway needs to be present in both datasets (i.e. conditions) to work", 
  file = log_file
)

for (pathway in shared_pathways) {

  pathways.show <- c(pathway)
  print(paste0("Visualizing pathway: ", pathway))

  # loop 


  tryCatch({
    
    # open device 
    png(
      filename = file.path(
        chord_comp_signalling_changes_de_dir, 
        paste0("chord-comp-",pathway,".png")
      ),
      width = 20, height = 10, res = 300, units = "in"
    )
    # plot
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(cellchat_list)) {
      
      netVisual_aggregate(
        cellchat_list[[i]], 
        signaling = pathways.show, 
        layout = "chord",
        signaling.name = paste(pathways.show, names(cellchat_list)[i])
      )
    }
    # close device
    dev.off()
    
  }, error = function(e) {
    
    msg <- paste(
      "Failed to visualize pathway:", 
      pathway, 
      "in dataset:", 
      names(cellchat_list)[i]
    )
    
    message(msg)
    

    # write to the log file 
    write(
      msg,
      file = log_file,
      append = TRUE
    )
  })


}

#######################################################################
# PART 5: COMPARE THE SIGNALLING GENE EXPRESSION DISTRIBUTION BETWEEN #
# DIFFERENT DATASETS                                                  #
#######################################################################

# set the dataset factor for the merged object, which is required for the
cellchat_merged@meta$datasets = factor(
  cellchat_merged@meta$condition, 
  levels = c("CTRL", "CF")
) 

# DIR for the violin plots for dysregulated pathways
violin_comp_signalling_changes_de_dir <- file.path(
  OUT_DIR, "violin-comp-signaling-changes-de"
)
dir.create(violin_comp_signalling_changes_de_dir, showWarnings = FALSE)

for (pathway in shared_pathways) {
  png(
    filename = file.path(
      violin_comp_signalling_changes_de_dir, 
      paste0("violin-", pathway, ".png")
    ),
    width = 10, height = 4, res = 300, units = "in"
  )
  print(
    plotGeneExpression(
    cellchat_merged, 
    signaling = pathway, 
    split.by = "datasets", 
    colors.ggplot = T, 
    type = "violin"
    )
  )
  dev.off()
}

# save the object list 
saveRDS(
  cellchat_list, 
  file.path(OUT_DIR, "cellchat_list.rds")
)

# save the merged object for future use
saveRDS(
  cellchat_merged, 
  file.path(OUT_DIR, "cellchat_merged.rds")
)

# save also the de object list and merged object for future use
saveRDS(
  cellchat_list_de, 
  file.path(OUT_DIR, "cellchat_list_de.rds")
)

saveRDS(
  cellchat_merged_de, 
  file.path(OUT_DIR, "cellchat_merged_de.rds")
)