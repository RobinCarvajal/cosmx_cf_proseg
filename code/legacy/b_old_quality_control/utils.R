
grub_test <- function(meta, col, type = 10, opposite = FALSE, two.sided = TRUE) {
  if (!requireNamespace("outliers", quietly = TRUE)) {
    stop("Please install the 'outliers' package: install.packages('outliers')")
  }

  if (!col %in% colnames(meta)) {
    stop(paste("Column", col, "not found in metadata"))
  }

  values <- meta[[col]]

  if (!is.numeric(values)) {
    stop("Selected column must be numeric.")
  }

  if (length(values) < 3) {
    stop("Grubbs' test requires at least 3 numeric values.")
  }

  result <- outliers::grubbs.test(values, type = type,
                                  opposite = opposite, two.sided = two.sided)
  return(result)
}



complexity_test <- function(meta, ratio = 1) {

  counts <- meta[["nCount_RNA"]]
  genes <- meta[["nFeature_RNA"]]

  meta$qc_complexity <- ifelse((counts / genes) >= ratio, "pass", "fail")

  return(meta)

}

negprobe_test <- function(meta, cutoff = 0.1) {

  prop_negprobes <- meta[["propNegative"]]
  meta$qc_prop_negprobes <- ifelse(prop_negprobes >= cutoff,
                                   "fail", "pass")

  return(meta)

}

count_test <- function(meta, cutoff = 20) {

  counts <- meta[["nCount_RNA"]]

  meta$qc_count <- ifelse(counts >= cutoff, "pass", "fail")

  return(meta)
}

area_test <- function(){
  next
}



cell_qc <- function(meta, complexity_ratio = 1, prop_negprobes_cutoff = 0.1,
                    count_cutoff = 20) {

  # complexity test
  meta <- complexity_test(meta = meta, ratio = complexity_ratio)

  # area test
  meta <- negprobe_test(meta = meta, cutoff = prop_negprobes_cutoff)

  # counts test
  meta <- count_test(meta = meta, cutoff = count_cutoff)

  


}