# %% Load libraries
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(pryr) # to check ram usage mem_used()

# %% Paths

# project dir
project_dir <- "/mnt/data/project0062/cosmx_gray"
# dirs
data_dir <- file.path(project_dir, "data")
diag_dir <- file.path(project_dir, "results", "diag")

# create directories to save the plots
rna_dir <- file.path(diag_dir, "rna")
# create
dir.create(rna_dir, recursive = TRUE)

# create directories to save the plots
negprobes_dir <- file.path(diag_dir, "negprobes")
# create
dir.create(negprobes_dir, recursive = TRUE)

# to save qc cols plots
qc_dir <- file.path(diag_dir, "qc")
# create
dir.create(qc_dir)

# files
rds_path <- file.path(data_dir, "combined/rds/Combined_v0.RDS")
samples_meta_path <- file.path(data_dir, "combined/meta/Combined_v1.csv")
new_rds_path <- file.path(data_dir, "combined/h5ad/Combined_v1.RDS")

# %% read obj
obj <- readRDS(rds_path)

# %% adding sample meta
samples_meta <- read.csv(samples_meta_path, row.names = 1)

# add meta back to obj
obj@meta.data <- samples_meta

# %% Remove cells where "not assigned" in sample
obj <- subset(obj, subset = sample != "not_assigned")

################
#### FOV QC ####
################

## source the necessary functions:
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

## load barcodes:
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))

barcodemap <- allbarcodes$Hs_6k

####### erase from here

# get the counts
test_sample <- "slide2"

smpl <- subset(obj, subset = slide_name == test_sample)
cts <- Matrix::t(smpl@assays$RNA@counts)
xy <- smpl@meta.data[c("x_slide_mm", "y_slide_mm")]
fov <- as.vector(smpl@meta.data["fov"])$fov
metadata <- smpl@meta.data

# parameters
max_prop_loss <- 0.6
max_totalcounts_loss <- 0.6

Sys.time()
# main function
res <- runFOVQC(
  counts = cts,
  xy = xy,
  fov = fov,
  tissue = paste0(metadata$sample, "_"),
  barcodemap = barcodemap,
  max_prop_loss = max_prop_loss,
  max_totalcounts_loss = max_totalcounts_loss
)
Sys.time()

# which FOVs were flagged
res$flaggedfovs
mapFlaggedFOVs(res)

head(res$flagged_fov_x_gene)
head(res$flagged_fov_x_gene[, "gene"])
res$flaggedfovs_fortotalcounts
res$flaggedfovs_forbias
# count how many genes were impacted in one or more flagged FOVs:
length(unique(res$flagged_fov_x_gene[, "gene"]))
FOVSignalLossSpatialPlot(res, shownames = TRUE)

library(pheatmap)
FOVEffectsHeatmap(res)

par(mfrow = c(4, 4))
FOVEffectsSpatialPlots(
  res = res,
  outdir = NULL,
  bits = "flagged_reportercycles"
)

par(mfrow = c(1, 1))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "reportercycle11R") 


flag <- is.element(metadata$fov, res$flaggedfovs)
table(flag)

# %% spliting the object by sample
obj$sample_name <- obj$sample
sample_names <- unique(obj$sample_name)
obj$sample <- NULL

sample_objs <- list()
# subset the obj for each sample and store in a list
for (sample in sample_names) {

  obj_sub <- subset(obj, subset = sample_name == sample)
  sample_objs[[sample]] <- obj_sub

}

# run form here ##################

init_time <- Sys.time()
# run FOV QC for each sample
results <- list()
for (sample in sample_names) {

  # info message
  print(paste("Processing sample:", sample))

  # assigning sample object
  smpl <- sample_objs[[sample]]

  # split the object in parts to be analysed
  cts <- Matrix::t(smpl@assays$RNA@counts)
  xy <- smpl@meta.data[c("x_slide_mm", "y_slide_mm")]
  metadata <- smpl@meta.data
  fov <- metadata$fov

  # parameters
  max_prop_loss <- 0.5
  max_totalcounts_loss <- 0.5

  # main function
  res <- runFOVQC(
    counts = cts,
    xy = xy,
    fov = fov,
    barcodemap = barcodemap,
    max_prop_loss = max_prop_loss,
    max_totalcounts_loss = max_totalcounts_loss
  )

  # append to results
  results[[sample]] <- res
}

final_time <- Sys.time()

total_seconds <- as.numeric(final_time - init_time, units = "secs")
cat(sprintf("FOV QC completed in %d minutes and %d seconds.\n",
            total_seconds %/% 60, round(total_seconds %% 60)))




res <- results$sample1

mapFlaggedFOVs(res)

# which FOVs were flagged
res$flaggedfovs
mapFlaggedFOVs(res)

head(res$flagged_fov_x_gene)
head(res$flagged_fov_x_gene[, "gene"])
res$flaggedfovs_fortotalcounts
res$flaggedfovs_forbias
# count how many genes were impacted in one or more flagged FOVs:
length(unique(res$flagged_fov_x_gene[, "gene"]))
FOVSignalLossSpatialPlot(res, shownames = TRUE)

library(pheatmap)
FOVEffectsHeatmap(res)

par(mfrow = c(4, 4))
FOVEffectsSpatialPlots(
  res = res,
  outdir = NULL,
  bits = "flagged_reportercycles"
)

par(mfrow = c(1, 1))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "reportercycle11R") 
