#!/usr/bin/env bash

em_table="/Volumes/robin_work/cosmx_gray/results/comb/ct_de/T/pb_T.csv"
ss_table="/Volumes/robin_work/cosmx_gray/data/sample_sheet/sample_sheet.tsv"
design_string='~ run + condition'
de_out_path="/Volumes/robin_work/cosmx_gray/results/comb/ct_de/T/de_T.tsv"
ne_out_path="/Volumes/robin_work/cosmx_gray/results/comb/ct_de/T/ne_T.tsv"

Rscript DESeq2.R --em $em_table --ss $ss_table --design "$design_string" \
  --de_out $de_out_path --ne_out $ne_out_path