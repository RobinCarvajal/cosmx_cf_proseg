# loading libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))

em_path = "/mnt/data/project0062/cosmx_gray/results/comb/ct_de/secretory-ciliated/pb_secretory-ciliated.csv"
ss_path = "/mnt/data/project0062/cosmx_gray/data/sample_sheet/sample_sheet.tsv"


# Load files
em <- read.csv(em_path, row.names = 1)
ss <- read.csv(ss_path, sep = '\t')
design_formula = ~ run + condition

# EXTEMELY IMPORTANT
# rearrange ss samples order based on em columns order
ss <- ss[match(colnames(em), ss$sample_name), , drop = FALSE]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = em, 
  colData = ss, 
  design = design_formula
)

# Filter lowly expressed genes
dds <- dds[rowMeans(counts(dds)) >= 1, ] # Keep genes with mean raw counts >= 1 

# setting  reference condition
dds$condition <- relevel(dds$condition, ref = "CTRL")

## Run DESEQ2
dds <- DESeq(dds)

#Checking coeficients
resultsNames(dds)

## results
res <- lfcShrink(dds,
                 coef = 4, # only use with 
                 type = "apeglm") # ashr or apeglm

# res_df <- searchlight_prep(res)

############ PCA BATCH CORRECTION

library(limma)

vsd <- vst(dds, blind = FALSE)  # uses the design (good for plots)

####### No Batch Correction ##########

p <- plotPCA(vsd, intgroup = c("condition", "run"), returnData = TRUE)

percentVar <- round(100 * attr(p, "percentVar"))

library(ggplot2)
library(ggrepel)
ggplot(p, aes(PC1, PC2, color = condition, shape = run)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = name), show.legend = FALSE, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%"))


####### With Batch Correction #########

mat <- assay(vsd)

mat_bc <- removeBatchEffect(
  em,
  batch = ss$run,
  design = model.matrix(~ condition, data = ss)
)

library(ggplot2)

# PCA on samples (so transpose: samples x genes)
pca <- prcomp(t(mat_bc), center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  batch = ss$run,
  condition = ss$condition
)

pca_df$condition <- factor(pca_df$condition, levels = c("CTRL", "CF"))

# % variance explained for axis labels
ve <- (pca$sdev^2) / sum(pca$sdev^2) * 100

ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_label_repel(aes(label = sample), show.legend = FALSE, size = 3) +
  xlab(sprintf("PC1 (%.1f%%)", ve[1])) +
  ylab(sprintf("PC2 (%.1f%%)", ve[2])) +
  scale_color_manual(values = c("CTRL" = "blue", "CF" = "red"))


########### save the batch corrected table ###########

design_terms <- attr(terms(design_formula), "term.labels")  # c("run","condition") etc.
batch_vec <- ss[[design_terms[1]]]
condition_formula <- as.formula(paste('~', design_terms[2]))

# setting  reference condition
ss$condition <- relevel(factor(ss$condition), ref = "CTRL")

vsd <- vst(dds, blind = FALSE)  # uses the design (good for plots)
mat <- assay(vsd)

bc_mat <- removeBatchEffect(
  em,
  batch = batch_vec,
  design = model.matrix(condition_formula, data = ss)
) 

bc_table <- data.frame(bc_mat)
bc_table$ID <- row.names(bc_table)
bc_table <- bc_table[, c("ID", setdiff(names(bc_table), "ID"))]

write.table(bc_table,
            file = args$bc_out,
            row.names = FALSE,
            quote = FALSE, 
            sep = "\t")


