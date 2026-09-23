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

# %% Defining samples
samples <- unique(obj@meta.data$sample)

# %% Defining rna parameters
rna_params <- grep(
  "rna", names(obj@meta.data),
  value = TRUE, ignore.case = TRUE
)

# %% [DIAG] Vln plots of rna columns

# save plots in a loop

for (rna_param in rna_params){
  vln_param <- VlnPlot(
    obj, features = rna_param,
    group.by = "sample",
    pt.size = 0
  )
  # plot path
  vln_plt_path <- file.path(
    rna_dir, paste0("vln_", rna_param, ".png")
  )
  # saving plot
  png(vln_plt_path, width = 10, height = 10, units = "in", res = 600)
  print(vln_param)
  dev.off()
}

# %% Defining nreprobes parameters
negprobes_params <- grep(
  "neg", names(obj@meta.data),
  value = TRUE, ignore.case = TRUE
)

# %% [DIAG] Vln plots of negprobes columns

# save plots in a loop
for (np_param in negprobes_params){
  plot <- VlnPlot(
    obj, features = np_param,
    group.by = "sample",
    pt.size = 0
  )
  # plot path
  plot_path <- file.path(
    negprobes_dir, paste0("vln_", np_param, ".png")
  )
  # saving plot
  png(plot_path, width = 10, height = 10, units = "in", res = 600)
  print(plot)
  dev.off()
}

# %% defining qc parameters

# get qc columns
qc_params <- grep("qc", names(samples_meta), value = TRUE, ignore.case = TRUE)
# remove qcCellsFlagged, its just the opposite of qcCellsPassed
qc_params <- qc_params[qc_params != "qcCellsFlagged"]

# %% Metrics in Barplots

# Loop through each QC flag
for (qc_param in qc_params) {

  # Create grouped summary: QC status by sample
  qc_table_df <- obj@meta.data |>
    count(sample, !!sym(qc_param)) |>
    rename(Status = !!sym(qc_param), Count = n)

  # Convert logical to character for clean facet labels
  qc_table_df$Status <- as.character(qc_table_df$Status)

  # Plot
  plot <- ggplot(qc_table_df, aes(x = sample, y = Count, fill = Status)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("QC Flag:", qc_param),
         x = "Sample", y = "Cell Count", fill = "QC Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(plot)

  # plot path
  plot_path <- file.path(
    qc_dir, qc_param, paste0("bp_", qc_param, ".png")
  )
  print(plot_path)
  # saving plot
  png(plot_path, width = 10, height = 10, units = "in", res = 600)
  print(plot)
  dev.off()

  # Remove this `break` if you want to plot all QC flags
  # break
}

# %% Metrics in Spatial Plots

# create directories to save the plots
for (qc_param in qc_params){
  qc_param_dir <- file.path(qc_dir, qc_param)
  # create
  dir.create(qc_param_dir, recursive = TRUE)
}

# save plots in a loop
for (slide in names(obj@images)){
  for (qc_param in qc_params){
    plot <- ImageDimPlot(
      obj, fov = slide,
      group.by = qc_param,
      size = 0.1,
      flip_xy = FALSE
    )
    # plot path
    sp_plt_path <- file.path(
      qc_dir, qc_param, paste0("sp_", slide, "_", qc_param, ".png")
    )
    # saving plot
    png(sp_plt_path, width = 10, height = 10, units = "in", res = 600)
    print(plot)
    dev.off()
  }
}


####################################################
#################### RANDOM STUFF ##################
####################################################

# change assay to negprobes
DefaultAssay(obj) <- "negprobes"

VlnPlot(obj, features = "nCount_negprobes", group.by = "sample")
VlnPlot(obj, features = "nFeature_negprobes", group.by = "sample")
# Other variables to check


sample_meta$nCount_falsecode


grep("qc", names(samples_meta), value = TRUE, ignore.case = TRUE)

grep("neg", names(samples_meta), value = TRUE, ignore.case = TRUE)
grep("prop", names(samples_meta), value = TRUE, ignore.case = TRUE)


DefaultAssay(obj) <- "falsecode"

VlnPlot(obj, features = "nCount_falsecode", group.by = "sample")
VlnPlot(obj, features = "nFeature_falsecode", group.by = "sample")


# %% check how many cells pass qc if nCount_RNA = 20 

meta <- obj@meta.data["nCount_RNA"]
min_counts <- 100
meta <- meta[meta$nCount_RNA >= min_counts,]
length(meta)

# summary
# 20 - 865227
# 50 - 767111
# 100 - 629173

obj_test <- obj

min_counts <- 50
obj_test@meta.data["nCount_RNA_20"] <-obj_test@meta.data["nCount_RNA"] >= min_counts # nolint

ImageDimPlot(
  obj_test,
  fov = "slide2",
  group.by = "nCount_RNA_20",
  flip_xy = FALSE
)
###



min_features <- 50
obj_test@meta.data["nFeature_RNA_50"] <- obj_test@meta.data["nFeature_RNA"] >= min_features # nolint

ImageDimPlot(
  obj_test,
  fov = "slide2",
  group.by = "nFeature_RNA_50",
  flip_xy = FALSE
)


grep("complex", names(samples_meta), value = TRUE, ignore.case = TRUE)

# kmeans to identify samples 
meta_slide1 <- obj@meta.data[obj@meta.data$slide_name == "slide1",]
coords <- meta_slide1[c("x_slide_mm", "y_slide_mm")]

# Run k-means clustering (e.g., with 3 clusters)
set.seed(123)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 2)

# Add cluster labels back to Seurat metadata
meta_slide1$spatial_kmeans <- factor(kmeans_result$cluster)

# Optional: visualize clusters on the spatial coordinates
library(ggplot2)
ggplot(meta_slide1, aes(x = x_slide_mm, y = y_slide_mm, color = spatial_kmeans)) +
  geom_point(size = 0.5) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "K-means Clustering of Spatial Coordinates")


  # kmeans to identify samples 
meta_slide2 <- obj@meta.data[obj@meta.data$slide_name == "slide2",]
coords <- meta_slide2[c("x_slide_mm", "y_slide_mm")]

# Run k-means clustering (e.g., with 3 clusters)
set.seed(1234)  # for reproducibility
kmeans_result <- kmeans(coords, centers = 3, iter.max = 20)

# Add cluster labels back to Seurat metadata
meta_slide2$spatial_kmeans <- factor(kmeans_result$cluster)

# Optional: visualize clusters on the spatial coordinates
library(ggplot2)
ggplot(meta_slide2, aes(x = x_slide_mm, y = y_slide_mm, color = spatial_kmeans)) +
  geom_point(size = 0.5) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "K-means Clustering of Spatial Coordinates")


###################
#### QUANTILES ####
###################

# meta for each slide
meta_slide1 <- obj@meta.data[obj@meta.data$slide_name == "slide1", ]
meta_slide2 <- obj@meta.data[obj@meta.data$slide_name == "slide2", ]

# q values
q_01_09 <- seq(0.1, 0.90, by = 0.1)
q_09_099 <- seq(0.9, 0.99, by = 0.01)

quantile(meta_slide2$unassignedTranscripts, probs = q_01_09)
sum(meta_slide2$nCount_RNA)
sum(meta_slide1$nCount_RNA)
sum(meta_slide2$nCount_RNA)/sum(meta_slide1$nCount_RNA)

quantile(meta_slide2$nCount_RNA, probs = q_01_09)
quantile(meta_slide2$median_negprobes, probs = q_01_09)

q_values <- seq(0.9, 0.99, by = 0.01)
quantile(meta_slide2$nCount_RNA, probs = q_values)

max(meta_slide2$nCount_RNA)
head(sort(meta_slide2$nCount_RNA, decreasing = TRUE), 5000)


samples_meta %>% glimpse()

hist(samples_meta$nCount_RNA, breaks = 50)

meta <- samples_meta[samples_meta$sample == "sample4", ]
hist(meta$nCount_RNA, breaks = 50)
plot(density(meta$nCount_RNA))
meta$log_nCount_RNA <- log1p(meta$nCount_RNA)  # log1p = log(x + 1)
hist(meta$log_nCount_RNA, breaks = 30)
plot(density(meta$log_nCount_RNA))

meta$log_nFeature_RNA <- log1p(meta$nFeature_RNA)  # log1p = log(x + 1)
hist(meta$log_nFeature_RNA, breaks = 30)
plot(density(meta$log_nFeature_RNA))

hist(meta$propNegative, breaks = 100)

quantile(meta$Area.um2, probs = c(0.01, 0.99))


percentiles <- quantile(meta$Area.um2,
                        probs = seq(0.02, 0.99, by = 0.01), na.rm = TRUE)
plot(density(percentiles))

result <- outliers::grubbs.test(percentiles, type = 11,
                                opposite = FALSE, two.sided = TRUE)
result
result$alternative
hist(meta$Area.um2, breaks = 100)

library(ggplot2)

# Create a simple data frame
df <- data.frame(value = meta$Area.um2)

# Add a dummy factor column for grouping
df$group <- "Area.um2"

# Plot
ggplot(df, aes(x = group, y = value)) +
  geom_violin(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3, color = "black") +
  labs(x = "", y = "Area (µm²)", title = "Violin plot with jitter") +
  theme_minimal()
# remove the top percentile

# Extract the alternative hypothesis string
alt_text <- result$alternative

# Use regex to extract the number from the string
outlier_value <- as.numeric(gsub(".*value\\s+(\\S+)\\s+is.*", "\\1", alt_text))


annotate_percentile <- function(meta, col = "nCount_RNA", 
                                probs = seq(0.01, 0.99, by = 0.01),
                                label_col = "percentile_bin") {
  if (!is.numeric(meta[[col]])) {
    stop("Selected column must be numeric.")
  }

  # Compute percentiles
  cuts <- quantile(meta[[col]], probs = probs, na.rm = TRUE)

  # Assign bin labels (e.g. "1%", "2%", ..., "99%", "100%")
  bin_labels <- c(names(cuts), "100%")

  # Cut values into bins (include Inf at the top)
  meta[[label_col]] <- cut(meta[[col]],
                           breaks = c(-Inf, cuts, Inf),
                           labels = bin_labels,
                           include.lowest = TRUE,
                           right = TRUE)

  return(meta)
}


perc_meta <- annotate_percentile(meta = samples_meta, col = "Area.um2")
obj@meta.data <- perc_meta
test <- subset(obj, subset = percentile_bin == "99%")
ImageDimPlot(test, fov = "slide2", 
             #group.by = "percentile_bin",
             cols = "red") + NoLegend()

test <- obj@meta.data[obj$slide_name == "slide2", ]

plot(test_tab$nCount_RNA, test_tab$Area.um2)
plot(test_tab$complexity, test_tab$Area.um2)
plot(test$unassignedTranscripts, test$slide_ID_numeric)
plot(test$nCount_RNA, test$slide_ID_numeric)
plot(test$Circularity, test$nCount_RNA)
plot(test$Circularity, test$slide_ID_numeric)
plot(test$Mean.Membrane, test$nCount_RNA)
plot(test$Mean.DAPI, test$nCount_RNA)

plot(test$nCount_RNA, test$Area.um2)
plot(test$nFeature_RNA, test$Area.um2)

plot(test$Area.um2, test$nFeature_RNA)

s1_meta <- obj@meta.data[obj@meta.data$sample == "sample4", ]
s1_meta <- annotate_percentile(meta = s1_meta, col = "Area.um2",
                               probs = seq(0.01, 0.99, by = 0.01))
library(ggplot2)
library(viridis)
ggplot(s1_meta, aes(x = x_slide_mm, y = y_slide_mm, colour = percentile_bin)) +
  geom_point() +
  scale_colour_viridis_d()


s1_meta <- obj@meta.data[obj@meta.data$sample == "sample3", ]
threshold <- quantile(s1_meta$nCount_RNA, 0.99, na.rm = TRUE)
s1_meta$top1percent <- s1_meta$nCount_RNA >= threshold
library(ggplot2)
library(viridis)
ggplot(s1_meta, aes(x = x_slide_mm, y = y_slide_mm, colour = top1percent)) +
  geom_point(size = 1) +
  scale_colour_manual(values = c("FALSE" = "grey80", "TRUE" = viridis(1))) +
  theme_minimal()


values <- s1_meta$Area
result <- outliers::grubbs.test(values, type = 10,
                                  opposite = FALSE, 
                                  two.sided = TRUE)

result


#### shit test 

library(outliers)

iterative_grubbs <- function(values, alpha = 0.05) {
  values <- na.omit(values)
  outliers <- c()
  
  repeat {
    test <- grubbs.test(values, type = 10, two.sided = TRUE)
    
    if (test$p.value < alpha) {
      # Find the most extreme value (Grubbs tests for max deviance from mean)
      deviation <- abs(values - mean(values))
      idx_to_remove <- which.max(deviation)
      outliers <- c(outliers, values[idx_to_remove])
      values <- values[-idx_to_remove]
    } else {
      break
    }
  }
  
  return(list(cleaned = values, outliers = outliers, final_test = test))
}

# Example usage
result <- iterative_grubbs(s1_meta$Area.um2)
result$cleaned        # values with outliers removed
result$outliers       # values detected as outliers
result$final_test     # last non-significant Grubbs test


hist(s1_meta$Area.um2, main = "Histogram of Area.um2", xlab = "Area.um2", breaks = 50)
abline(v = min(result$outliers), col = "red", lwd = 2)
abline(v = 250 , col = "blue", lwd = 2)
abline(v = 25 , col = "green", lwd = 2)



percentiles <- quantile(s1_meta$Area.um2,
                        probs = seq(0.01, 0.99, by = 0.01), na.rm = TRUE)

small_cells <- s1_meta[s1_meta$Area.um2 <= 25, ]
plot(small_cells$nCount_RNA, small_cells$Area.um2)
plot(small_cells$nFeature_RNA, small_cells$Area.um2)

plot(density(small_cells$nCount_RNA))

library(dplyr)
small_cells %>% select(nCount_RNA, Area.um2) 


lowqc_cells <- s1_meta[s1_meta$nCount_RNA <= 50, ]
plot(lowqc_cells$Area.um2, lowqc_cells$nCount_RNA)
plot(lowqc_cells$Area.um2, lowqc_cells$nFeature_RNA)

plot(density(small_cells$nCount_RNA))

library(dplyr)
lowqc_cells %>% select(nCount_RNA, Area.um2) 



high_exp_cells <- samples_meta[samples_meta$nCount_RNA >= 500, ]
plot(high_exp_cells$Area.um2, high_exp_cells$nCount_RNA)
plot(high_exp_cells$Area.um2, high_exp_cells$nFeature_RNA)

plot(density(high_exp_cells$nCount_RNA))

library(dplyr)
high_exp_cells %>% select(nCount_RNA, Area.um2) 


### get cells that thouch the borders of the fov

# get fov border coordinates
for 


quantile(samples_meta$Solidity,
         probs = seq(0.01, 0.99, by = 0.01),
         na.rm = TRUE)

quantile(samples_meta$Circularity,
         probs = seq(0.01, 0.99, by = 0.01),
         na.rm = TRUE)


samples_meta$Solidity_corrected <- 1 / samples_meta$Solidity
quantile(samples_meta$Solidity_corrected,
         probs = seq(0.01, 0.99, by = 0.01),
         na.rm = TRUE)
glimpse(samples_meta)

quantile(samples_meta$nFeature_RNA,
         probs = seq(0.0, 1, by = 0.01),
         na.rm = TRUE)

quantile(samples_meta$nCount_RNA,
         probs = seq(0.0, 1, by = 0.01),
         na.rm = TRUE)
