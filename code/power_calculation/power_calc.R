X <- read.csv("/Volumes/robin_work/cosmx_gray/results/comb/power/ciliated/vst_table.tsv", sep='\t', row.names=1)
meta <- read.csv("/Volumes/robin_work/cosmx_gray/data/sample_sheet/sample_sheet.tsv", sep='\t')
de_table <- read.csv("/Volumes/robin_work/cosmx_gray/results/comb/power/ciliated/de_table.tsv", sep='\t', row.names=1)

# rearrange ss samples order based on em columns order
meta <- meta[match(colnames(X), meta$donor), , drop = FALSE] # can be sample_name as well 

cohens_d <- function(xA, xC){
  nA <- length(xA); nC <- length(xC)
  sP <- sqrt(((nA-1)*var(xA) + (nC-1)*var(xC)) / (nA+nC-2))
  (mean(xA) - mean(xC)) / sP
}

condA <- meta$condition == "CF"   # change if your label is CF
condC <- meta$condition == "CTRL"   # change if your label is CTRL

d_vec <- apply(X, 1, \(g) cohens_d(g[condA], g[condC]))

library(pwr)

# de_table must have rownames = genes, and a column "padj"
top_genes <- rownames(de_table[order(de_table$P.adj), ])[1:2000] # 2000 top genes

d_med <- median(abs(d_vec[top_genes]), na.rm=TRUE)

pwr.t.test(d = d_med, power = 0.80, sig.level = 0.05, type = "two.sample")$n


d_cons <- quantile(abs(d_vec[top_genes]), 0.25, na.rm=TRUE)
pwr.t.test(d = d_cons, power=0.8, sig.level=0.05, type="two.sample")$n



