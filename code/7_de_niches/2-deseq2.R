#!/usr/bin/env Rscript

# loading libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(limma))

# searchlight prep function
searchlight_prep <- function(res){
  
  # converting to a df
  res_df <- data.frame(res)
  
  # order by padj
  res_df <- res_df[order(res_df$padj),]
  # gene(ID) column
  res_df$ID <- rownames(res_df)
  # Move 'ID' column to the first position
  res_df <- res_df[, c("ID", setdiff(names(res_df), "ID"))]
  
  # subset columns
  res_df <- select(res_df, ID, log2FoldChange, pvalue, padj)
  # rename the columns
  colnames(res_df) <- c("ID", "log2Fold", "P", "P.adj")
  
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
dds$condition <- relevel(dds$condition, ref = ref)

## Run DESEQ2
dds <- DESeq(dds)

#Checking coeficients
resultsNames(dds)


## results
res <- lfcShrink(dds,
                 coef = coef_string, # only use with 
                 type = "apeglm") # ashr or apeglm

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
group_formula <- as.formula(paste('~', group))

# setting  reference condition
ss$group <- relevel(factor(ss[[group]]), ref = ref)

# had to use a trick here for the very small niches
m <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
print(sum(m > 5))  # how many genes qualify?

vsd <- vst(dds, blind = FALSE, nsub=150)  # uses the design (good for plots)
mat <- assay(vsd)

bc_mat <- removeBatchEffect(
  em,
  batch = batch_vec,
  design = model.matrix(group_formula, data = ss)
) 

bc_table <- data.frame(bc_mat)
bc_table$ID <- row.names(bc_table)
bc_table <- bc_table[, c("ID", setdiff(names(bc_table), "ID"))]

write.table(bc_table,
            file = args$bc_out,
            row.names = FALSE,
            quote = FALSE, 
            sep = "\t")