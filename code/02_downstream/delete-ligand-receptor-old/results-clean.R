# Load libraries and scripts ---------------------------------------------------

library(Seurat)
library(schard)
#library(SpatialCellChat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
source("code/8_lr_analysis/cellchat-functions.R")

# Set paths --------------------------------------------------------------------
MAIN_DIR <- "/mnt/autofs/data/userdata/project0062/cosmx_gray"
setwd(MAIN_DIR)

# Setting seed
SEED_VALUE <- 42
set.seed(SEED_VALUE)

# Set up parallel processing ---------------------------------------------------
future::plan("multisession", workers = 8) 
options(future.globals.maxSize = 4 * 1024^3)
# helps speed up LR identification

# Niche object path ------------------------------------------------------------

IN_OBJ_PATH <- file.path(
  MAIN_DIR, "results", "comb", "lr_manual-niches", "1", "CF.cellchat.rds"
)

OUT_DIR <- file.path(
  MAIN_DIR, "results", "comb", "lr_manual-niches", "1"
)

# load the niche object --------------------------------------------------------
cc_obj <- readRDS(IN_OBJ_PATH)

# ---------------------------------

interactions <- subsetCommunication(cc_obj)

groupSize <- as.numeric(table(cc_obj@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc_obj@net$count, vertex.weight = rowSums(cc_obj@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_obj@net$weight, vertex.weight = rowSums(cc_obj@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



basal_interactions <- interactions[
  interactions$source == "basal" | interactions$target == "basal",
]

# save interactions and basal interactions table -----------------------------

write.csv(interactions, file.path(OUT_DIR, "CF_interactions.csv"), row.names = FALSE)
write.csv(basal_interactions, file.path(OUT_DIR, "CF_basal_interactions.csv"), row.names = FALSE)

# -------------------------------------------

pathways.show <- unique(basal_interactions$pathway_name)

for (pathway in pathways.show) {

  png(
    filename = file.path(OUT_DIR, sprintf("circle_plot_%s.png", pathway)),
    width = 6,
    height = 6,
    units="in",
    res = 300
  )

  par(mfrow = c(1, 1), xpd = TRUE)

  netVisual_aggregate(
    cc_obj,
    signaling = pathway,
    layout = "circle"
  )

  dev.off()
}


#########

# Compute the network centrality scores
cc_obj <- netAnalysis_computeCentrality(cc_obj, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

for (pathway in pathways.show){
  print(pathway)

  png(
    filename = file.path(OUT_DIR, sprintf("signaling_role_network_%s.png", pathway)),
    width = 6,
    height = 3,
    units = "in",
    res = 300
  )

  par(mfrow=c(1,1))
  
  netAnalysis_signalingRole_network(
    cc_obj, 
    signaling = pathway, 
    width = 8, 
    height = 2.5, 
    font.size = 10
  )

  dev.off()
}
