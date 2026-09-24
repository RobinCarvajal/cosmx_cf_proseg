#!/usr/bin/env Rscript

# loading libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(dplyr))


# searchlight prep function
searchlight_prep <- function(res) {

  # Convert DESeq2 results to data frame
  res_df <- as.data.frame(res)

  # Add gene ID column
  res_df$ID <- rownames(res_df)

  # Move ID to the first column
  res_df <- res_df[, c("ID", setdiff(names(res_df), "ID"))]

  # Rename columns
  res_df <- res_df |>
    dplyr::rename(
      log2Fold = log2FoldChange,
      P = pvalue,
      P.adj = padj
    ) |>
    dplyr::arrange(P.adj)

  return(res_df)
}

parser <- ArgumentParser(description = "DESeq2 CLI")

parser$add_argument("--em",  required = TRUE, help = "Expression Matrix file")
parser$add_argument("--ss", required = TRUE, help = "Sample Sheet file")
parser$add_argument("--ref", required = TRUE, help = "Reference group e.g. CTRL")
parser$add_argument("--query", required = TRUE, help = "Query group e.g. TREATMENT")
parser$add_argument("--batch", required = TRUE, help = "Batch variable e.g. RUN")
parser$add_argument("--group", required = TRUE, help = "Group variable e.g. CONDITION")
parser$add_argument("--de_out", required = TRUE, help = "Path for the Diferential Expression table")
parser$add_argument("--ne_out", required = TRUE, help = "Path for the Nornalised Expression table")
parser$add_argument("--bc_out", required = TRUE, help = "Path for the Batch Corrected Expression table")

args <- parser$parse_args()

cat("em  =", args$em,  "\n")
cat("ss =", args$ss, "\n")

if (!file.exists(args$em)) stop("Input not found: ", args$em)
if (!file.exists(args$ss)) stop("Input not found: ", args$ss)

# Load files
em <- read.csv(args$em, row.names = 1)
ss <- read.csv(args$ss, sep = '\t')

## argument values
batch <- args$batch
group <- args$group
ref <- args$ref
query <- args$query

design_string <- paste0("~ ",batch," + ",group)
design_formula <- as.formula(design_string)

print(design_formula) # to console

coef_string <- paste0(group, "_", query,"_vs_", ref)

# EXTEMELY IMPORTANT
# rearrange ss samples/donor order based on em columns order
ss <- ss[match(colnames(em), ss$donor), , drop = FALSE]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = em, 
  colData = ss, 
  design = design_formula
)

# Filter lowly expressed genes
dds <- dds[rowMeans(counts(dds)) >= 1, ] # Keep genes with mean raw counts >= 1 

# setting  reference condition
dds$condition <- relevel(dds$condition, ref = ref)

## Run DESEQ2
dds <- DESeq(dds)

#Checking coeficients
resultsNames(dds)


## results
# res <- lfcShrink(dds,
#                  coef = coef_string, # only use with 
#                  type = "apeglm") # ashr or apeglm

res <- results(dds)

res_df <- searchlight_prep(res)

# save the ne_table
ne_table <- data.frame(counts(dds, normalized=TRUE))
ne_table$ID <- row.names(ne_table)
ne_table <- ne_table[, c("ID", setdiff(names(ne_table), "ID"))]
write.table(ne_table,
            file = args$ne_out,
            row.names = FALSE,
            quote = FALSE, 
            sep = "\t")

# save the de_table
write.table(res_df,
            file = args$de_out,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Get batch corrected expression matrix ----------------------------------------

batch_vec <- ss[[batch]]

# Set the reference condition
ss$group <- relevel(
  factor(ss[[group]]),
  ref = ref
)

# The design must use the newly created ss$group column
group_formula <- ~ group

# ---------------------------------------------------------
# Variance-stabilising transformation with automatic fallback
# ---------------------------------------------------------
normalized_counts <- counts(dds, normalized = TRUE)

n_eligible <- sum(
  rowMeans(normalized_counts, na.rm = TRUE) > 5
)

message(
  "[INFO] Genes with mean normalized count > 5: ",
  n_eligible
)

if (n_eligible >= 1000) {

  message("[INFO] Using fast vst()")

  vsd <- vst(
    dds,
    blind = FALSE,
    nsub = 1000
  )

} else {

  message(
    "[INFO] Fewer than 1000 eligible genes; ",
    "using varianceStabilizingTransformation()"
  )

  vsd <- varianceStabilizingTransformation(
    dds,
    blind = FALSE
  )
}

mat <- assay(vsd)

# ---------------------------------------------------------
# Remove batch effect from VST-transformed expression
# while preserving condition differences
# ---------------------------------------------------------
bc_mat <- removeBatchEffect(
  mat,
  batch = batch_vec,
  design = model.matrix(group_formula, data = ss)
)

# ---------------------------------------------------------
# Export
# ---------------------------------------------------------
bc_table <- data.frame(
  ID = rownames(bc_mat),
  bc_mat,
  check.names = FALSE
)

write.table(
  bc_table,
  file = args$bc_out,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)