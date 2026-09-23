# Set Working Dir ---------------------------------------------------------
MAIN_DIR <- "/mnt/autofs/data/userdata/project0062/cosmx_gray"
MAIN_DIR <- "/Volumes/robin_work/cosmx_gray"
setwd(MAIN_DIR)

# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
#library(SpatialCellChat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
source("code/8_lr_analysis/cellchat-functions.R")

# Set paths --------------------------------------------------------------------


# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Niche object path ------------------------------------------------------------
IN_DIR <- file.path(
  MAIN_DIR, "results", "comb", "lr_manual-niches", "1"
)

OUT_DIR <- file.path(
  MAIN_DIR, "results", "comb", "lr_manual-niches", "1", "comparison"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Define input CellChat objects ------------------------------------------------
# Adjust these names if needed
obj_paths <- c(
  CTRL = file.path(IN_DIR, "CTRL.cellchat.rds"),
  CF   = file.path(IN_DIR, "CF.cellchat.rds")
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
object_list <- lapply(obj_paths, readRDS)

# Ensure names are preserved
names(object_list) <- names(obj_paths)

for (item in names(object_list)) {
  # calculate centrality scores if not already done
  object_list[[item]] <- netAnalysis_computeCentrality(
    object_list[[item]],
    slot.name = "netP"
  )
}


# Merge CellChat objects -------------------------------------------------------
message("Merging CellChat objects...")
cellchat_merged <- mergeCellChat(
  object.list = object_list,
  add.names = names(object_list)
)

saveRDS(
  cellchat_merged,
  file = file.path(OUT_DIR, "cellchat_merged.rds")
)

# Basic comparison: number and strength of interactions ------------------------
# Number of interactions
png(
  filename = file.path(OUT_DIR, "compareInteractions_count.png"),
  width = 1800, height = 900, res = 150
)
compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
dev.off()

# Interaction strength
png(
  filename = file.path(OUT_DIR, "compareInteractions_weight.png"),
  width = 1800, height = 900, res = 150
)
compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)
dev.off()


# Differential interaction networks -------------------------------------------
# Red = increased in second dataset in group comparison
# Blue = decreased in second dataset in group comparison

png(
  filename = file.path(OUT_DIR, "diffInteraction_count.png"),
  width = 1600, height = 1600, res = 180
)
netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "count",
  label.edge = FALSE
)
dev.off()

png(
  filename = file.path(OUT_DIR, "diffInteraction_weight.png"),
  width = 1600, height = 1600, res = 180
)
netVisual_diffInteraction(
  cellchat_merged,
  weight.scale = TRUE,
  measure = "weight",
  label.edge = FALSE,
  #sources.use = c(16),
)
dev.off()

# Heatmaps of differential interactions ---------------------------------------
png(
  filename = file.path(OUT_DIR, "heatmap_diff_count.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_heatmap(cellchat_merged, measure = "count")
dev.off()

png(
  filename = file.path(OUT_DIR, "heatmap_diff_weight.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_heatmap(cellchat_merged, measure = "weight")
dev.off()

# Compare the number of interactions and interaction strength among different cell populations

# ((C) Circle plot showing the number of interactions or interaction
#  strength among different cell populations across multiple datasets

weight.max <- getMaxWeight(object_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object_list)) {
  netVisual_circle(object_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object_list)[i]))
}

# D

group.cellType <- c(
  "MYELOID",      # alveolar_macrophages
  "OTHER",        # ambiguous
  "EPITHELIAL",   # AT1
  "EPITHELIAL",   # AT2
  "LYMPHOID",     # B
  "EPITHELIAL",   # basal
  "EPITHELIAL",   # ciliated
  "ENDOTHELIAL",  # endothelial
  "STROMAL",      # fibroblasts
  "MYELOID",      # macrophages
  "MYELOID",      # mast
  "MYELOID",      # neutrophils
  "LYMPHOID",     # plasma
  "EPITHELIAL",   # secretory
  "EPITHELIAL",   # secretory-ciliated
  "STROMAL",      # smooth_muscle
  "OTHER",        # sputum
  "LYMPHOID"      # T
)

group.cellType <- factor(
  group.cellType,
  levels = c("EPITHELIAL", "MYELOID", "LYMPHOID", "ENDOTHELIAL", "STROMAL", "OTHER")
)

object_list_grouped <- lapply(object_list, function(x) {
  mergeInteractions(x, group.cellType)
})

cellchat_grouped <- mergeCellChat(
  object_list_grouped,
  add.names = names(object_list_grouped)
)

weight.max <- getMaxWeight(
  object_list_grouped,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "count", "count.merged")
)

png(
  file.path(OUT_DIR, "circle_grouped_number_interactions.png"),
  width = 2400, height = 1200, res = 180
)

par(mfrow = c(1, length(object_list_grouped)), xpd = TRUE)

for (i in seq_along(object_list_grouped)) {
  netVisual_circle(
    object_list_grouped[[i]]@net$count.merged,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = weight.max[3],
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", names(object_list_grouped)[i])
  )
}

dev.off()

weight.max.str <- getMaxWeight(
  object_list_grouped,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "weight", "weight.merged")
)

png(
  file.path(OUT_DIR, "circle_grouped_interaction_strength.png"),
  width = 2400, height = 1200, res = 180
)

par(mfrow = c(1, length(object.list_grouped)), xpd = TRUE)

for (i in seq_along(object_list_grouped)) {
  netVisual_circle(
    object_list_grouped[[i]]@net$weight.merged,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.weight.max = weight.max.str[3],
    edge.width.max = 12,
    title.name = paste0("Interaction strength - ", names(object_list_grouped)[i])
  )
}

dev.off()

# Comapare the major sources and targets in a 2D space -------------------------

num.link <- sapply(object_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object_list[[i]], title = names(object_list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


# (B) Identify the signaling changes of specific cell populations


celltypes <- levels(cellchat_merged@meta$ann_lvl_3)

for (celltype in celltypes) {

  # if fails next celltype
  tryCatch({
    gg <- netAnalysis_signalingChanges_scatter(
      cellchat_merged, idents.use = celltype, 
      #signaling.exclude = "MIF"
    )
    print(gg)
  }, error = function(e) {
    message("Failed for celltype: ", celltype)
  })
}

# Signalling groups

cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat_merged <- netClustering(cellchat_merged, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)




rankSimilarity(cellchat_merged, type = "functional")

gg1 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

#####


library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object_list[[i]]@netP$pathways, object_list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object_list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object_list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object_list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object_list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


netVisual_bubble(
  cellchat_merged, 
  #sources.use = 4, 
  #targets.use = c(5:11),  
  comparison = c(1, 2), 
  angle.x = 45
)

for (i in 1:length(celltypes)) {
  gg <- netVisual_bubble(
    cellchat_merged, 
    sources.use = i,
    comparison = c(1, 2),
    angle.x = 45,
    title.name = paste0("Bubble plot - ", celltypes[i])
  )
  print(gg)
}




gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

for (i in 1:length(celltypes)) {

    tryCatch({
    gg1 <- netVisual_bubble(
      cellchat_merged, 
      sources.use = i, 
      comparison = c(1, 2), 
      max.dataset = 2, 
      title.name = paste0("Increased signaling in CF - ", celltypes[i]), 
      angle.x = 45, 
      remove.isolate = TRUE
    )
    
    gg2 <- netVisual_bubble(
      cellchat_merged, 
      sources.use = i, 
      comparison = c(1, 2), 
      max.dataset = 1, 
      title.name = paste0("Decreased signaling in CF - ", celltypes[i]), 
      angle.x = 45, 
      remove.isolate = TRUE
    )
    
    print(gg1 + gg2)
  }, error = function(e) {
    message("Failed for celltype: ", celltypes[i])
  })
}


#### de analysis

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "CF"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged, 
  group.dataset = "datasets", 
  pos.dataset = pos.dataset, 
  features.name = features.name, 
  only.pos = FALSE, 
  thresh.pc = 0.1, 
  thresh.fc = 0.05,
  thresh.p = 0.05, 
  group.DE.combined = FALSE
) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(
  cellchat_merged, 
  features.name = features.name,
   variable.all = TRUE
)
# extract the ligand-receptor pairs with upregulated ligands in CF
net.up <- subsetCommunication(
  cellchat_merged, net = net, 
  datasets = "CF",
  ligand.logFC = 0.05, 
  receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in CTRL, i.e.,downregulated in CF
net.down <- subsetCommunication(
  cellchat_merged, net = net, 
  datasets = "CTRL",
  ligand.logFC = 0.05, 
  receptor.logFC = NULL
)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merged)

# A bubble plot

pairLR.use.up = net.up[, "interaction_name", drop = F]
sources.up <- unique(net.up$source)
sources.up.idx <- which(celltypes %in% sources.up)

for (source in sources.up.idx) {
  message("Processing source cell type: ", celltypes[source])

  net.up.celltype <- net.up[net.up$source == celltypes[source], ]
  pairLR.use.up <- net.up.celltype[, "interaction_name", drop = F]
  targets.up <- unique(net.up.celltype$target)
  targets.up.idx <- which(celltypes %in% targets.up)

  gg <- netVisual_bubble(
    cellchat_merged, 
    pairLR.use = pairLR.use.up, 
    sources.use = source, 
    targets.use = targets.up.idx,
    comparison = c(1, 2),  
    angle.x = 90, 
    remove.isolate = TRUE,
    title.name = paste0("Up-regulated signaling in ", names(object_list)[2], " - Source: ", celltypes[source])
  )
  print(gg)

}


gg1 <- netVisual_bubble(
  cellchat_merged, 
  pairLR.use = pairLR.use.up, 
  sources.use = sources.up.idx, 
  comparison = c(1, 2),  
  angle.x = 90, 
  remove.isolate = TRUE,
  title.name = paste0("Up-regulated signaling in ", names(object_list)[2])
)

pairLR.use.down = net.down[, "interaction_name", drop = F]
sources.down <- unique(net.down$source)
sources.down.idx <- which(celltypes %in% sources.down)
targets.down <- unique(net.down$target)
targets.down.idx <- which(celltypes %in% targets.down)
gg2 <- netVisual_bubble(
  cellchat_merged, 
  pairLR.use = pairLR.use.down, 
  sources.use = sources.down.idx, 
  targets.use = targets.down.idx,
  comparison = c(1, 2),  
  angle.x = 90, 
  remove.isolate = TRUE,
  title.name = paste0("Down-regulated signaling in ", names(object_list)[2])
)

print(gg1 + gg2)


##### net.up and net.down have the infromation of sender and receptor and lr pair
# it just need to be include din the function, which turns out is a general function for subsetting the communication network based on any criteria, which is very useful for downstream analysis.


# Chord diagram - looks horrible, doesnt give much infromation 


# this wordcloud is not very informative, 
# but it can be used to quickly check the most frequent ligands/receptors 
# in the up and down regulated interactions.

computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

## get up pathways 
pathways.up <- unique(net.up[,"pathway_name", drop = F])[,1]

# visually compare cell-cell communication - VERY INFROAMTIVE, at least visually
# not quantitatively

for (pathway in pathways.up) {

  pathways.show <- c(pathway)
  print(paste0("Visualizing pathway: ", pathway))

  tryCatch({
    weight.max <- getMaxWeight(
      object_list, 
      slot.name = c("netP"), 
      attribute = pathways.show
    ) # control the edge weights across different datasets
    
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object_list)) {
      netVisual_aggregate(
        object_list[[i]], 
        signaling = pathways.show, 
        layout = "circle", 
        edge.weight.max = weight.max[1], 
        edge.width.max = 10, 
        signaling.name = paste(pathways.show, names(object_list)[i]))
    }
  }, error = function(e) {
    message("Failed to visualize pathway: ", pathway)
  })
}

pathways.show <- c("CEACAM") # example pathways to show, adjust as needed
weight.max <- getMaxWeight(object_list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object_list)) {
  netVisual_aggregate(object_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object_list)[i]))
}

#3 this does not work unles the pathways is in both datasets 

pathways.show <- c("CEACAM") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object_list)) {
  ht[[i]] <- netVisual_heatmap(object_list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object_list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))



# Chord diagram
pathways.show <- c("CEACAM") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object_list)) {
  netVisual_aggregate(object_list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object_list)[i]))
}



cellchat_merged@meta$datasets = factor(cellchat_merged@meta$condition, levels = c("CTRL", "CF")) # set factor level
plotGeneExpression(cellchat_merged, signaling = "CEACAM", split.by = "datasets", colors.ggplot = T, type = "violin")

plotGeneExpression(cellchat_merged, signaling = "CLDN", split.by = "datasets", colors.ggplot = T, type = "violin")

## kinda works but its not very good, violin pltos never are anyway
# Ranking signaling pathways ---------------------------------------------------
png(
  filename = file.path(OUT_DIR, "rankNet_stacked.png"),
  width = 2000, height = 1400, res = 180
)
rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)
dev.off()

png(
  filename = file.path(OUT_DIR, "rankNet_unstacked.png"),
  width = 2000, height = 1400, res = 180
)
rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = FALSE,
  do.stat = TRUE
)
dev.off()

# Rank pathways specifically by information flow ------------------------------
# This can help separate shared vs condition-specific pathways
pdf(
  file = file.path(OUT_DIR, "rankNet_pairwise.pdf"),
  width = 12, height = 8
)
rankNet(
  cellchat_merged,
  mode = "comparison",
  comparison = c(1, 2),
  stacked = TRUE,
  do.stat = TRUE
)
dev.off()

# Compare outgoing / incoming signaling roles ---------------------------------
# Scatter plots for signaling changes at cell-group level
pdf(
  file = file.path(OUT_DIR, "signalingRole_scatter_all.pdf"),
  width = 10, height = 8
)
netAnalysis_signalingRole_scatter(cellchat_merged)
dev.off()

# Identify conserved and context-specific signaling ----------------------------
# Functional similarity requires same cell population composition across datasets
message("Computing network similarity and manifold embedding...")
cellchat_merged <- computeNetSimilarityPairwise(
  cellchat_merged,
  type = "functional"
)

cellchat_merged <- netEmbedding(
  cellchat_merged,
  type = "functional"
)

cellchat_merged <- netClustering(
  cellchat_merged,
  type = "functional"
)

cellchat_merged <- netAnalysis_computeCentrality(
  cellchat_merged,
  slot.name = "netP"
)

saveRDS(
  cellchat_merged,
  file = file.path(OUT_DIR, "cellchat_merged_similarity.rds")
)

png(
  filename = file.path(OUT_DIR, "netEmbedding_functional.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_embedding(
  cellchat_merged,
  type = "functional",
  label.size = 3.5
)
dev.off()

png(
  filename = file.path(OUT_DIR, "netEmbeddingZoomIn_functional.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_embeddingZoomIn(
  cellchat_merged,
  type = "functional",
  nCol = 2
)
dev.off()

# Structural similarity can also be useful if compositions differ slightly ----
message("Computing structural similarity and manifold embedding...")
cellchat_merged <- computeNetSimilarityPairwise(
  cellchat_merged,
  type = "structural"
)

cellchat_merged <- netEmbedding(
  cellchat_merged,
  type = "structural"
)

cellchat_merged <- netClustering(
  cellchat_merged,
  type = "structural"
)

saveRDS(
  cellchat_merged,
  file = file.path(OUT_DIR, "cellchat_merged_similarity_structural.rds")
)

png(
  filename = file.path(OUT_DIR, "netEmbedding_structural.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_embedding(
  cellchat_merged,
  type = "structural",
  label.size = 3.5
)
dev.off()

png(
  filename = file.path(OUT_DIR, "netEmbeddingZoomIn_structural.png"),
  width = 1800, height = 1400, res = 180
)
netVisual_embeddingZoomIn(
  cellchat_merged,
  type = "structural",
  nCol = 2
)
dev.off()

# Pathway-specific comparisons -------------------------------------------------
# Choose pathways present in merged object
pathways_all <- cellchat_merged@netP$pathways
writeLines(pathways_all, con = file.path(OUT_DIR, "pathways_detected.txt"))

# Example: inspect top pathways manually later
# Uncomment and edit if you want specific pathways
# pathways_use <- c("MIF", "CXCL", "TGFb")

# for (pathway in pathways_use) {
#   png(
#     filename = file.path(OUT_DIR, paste0("bubble_", pathway, ".png")),
#     width = 1800, height = 1400, res = 180
#   )
#   netVisual_bubble(
#     cellchat_merged,
#     signaling = pathway,
#     comparison = c(1, 2),
#     angle.x = 45
#   )
#   dev.off()
# }

# Save session info ------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUT_DIR, "sessionInfo.txt")
)

  message("Comparison analysis finished.")