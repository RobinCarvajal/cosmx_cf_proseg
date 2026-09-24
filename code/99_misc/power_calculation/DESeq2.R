#!/usr/bin/env Rscript

# loading libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))

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
parser$add_argument("--design", required = TRUE, help = "DESeq2 Design")
parser$add_argument("--de_out", required = TRUE, help = "Path for the DE table")
parser$add_argument("--ne_out", required = TRUE, help = "Path for the NE table")

args <- parser$parse_args()

cat("em  =", args$em,  "\n")
cat("ss =", args$ss, "\n")
cat("design =", args$ss, "\n")

if (!file.exists(args$em)) stop("Input not found: ", args$em)
if (!file.exists(args$ss)) stop("Input not found: ", args$ss)

# Load files
em <- read.csv(args$em, row.names = 1)
ss <- read.csv(args$ss, sep = '\t')
design_formula = as.formula(args$design)

# EXTEMELY IMPORTANT
# rearrange ss samples order based on em columns order
ss <- ss[match(colnames(em), ss$donor), , drop = FALSE] # can be sample_name as well 

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

res_df <- searchlight_prep(res)

# save the ne_table
# ne_table <- data.frame(counts(dds, normalized=TRUE))
# ne_table$ID <- row.names(ne_table)
# ne_table <- ne_table[, c("ID", setdiff(names(ne_table), "ID"))]

vsd <- vst(dds, blind=FALSE)
X <- assay(vsd)  # genes x donors
X <- data.frame(X)
X$ID <- row.names(X)
X <- X[, c("ID", setdiff(names(X), "ID"))]

write.table(X,
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
