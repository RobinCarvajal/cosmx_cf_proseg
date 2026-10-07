library(whereami)
library(qs)

wd <- dirname(thisfile())
setwd(wd)

## Map detailed cell-type labels to broader groups
label_to_group <- c(
  # AIRWAY EPITHELIUM
  BC      = "AIR_EPI",
  Club    = "AIR_EPI",
  MCC     = "AIR_EPI",
  Gob     = "AIR_EPI",
  Hillock = "AIR_EPI",
  Deu     = "AIR_EPI",

  # ALVEOLAR EPITHELIUM
  AT1 = "ALV_EPI",
  AT2 = "ALV_EPI",

  # LYMPHOID
  CD4_T     = "LYM",
  CD8_T     = "LYM",
  cLym      = "LYM",
  B         = "LYM",
  cB        = "LYM",
  PC_K_IgA  = "LYM",
  PC_K_IgG  = "LYM",
  PC_K_IgM  = "LYM",
  PC_L_IgA  = "LYM",
  PC_L_IgG  = "LYM",
  PC_K_IgAG = "LYM",

  # MYELOID
  Neut   = "MYE",
  Mast   = "MYE",
  Mac    = "MYE",
  AlvMac = "MYE",

  # STROMAL
  Fib    = "STROMA",
  actFib = "STROMA",
  SMC    = "STROMA",

  # ENDOTHELIAL
  EC    = "EC",
  lymEC = "EC",

  # AMB
  Amb = "AMB"

)

# all groups
all_groups <- unique(unname(label_to_group))

## Store both objects together
celltype_config <- list(
  label_to_group = label_to_group,
  all_groups = all_groups
)

## Save configuration
qsave(
  celltype_config,
  file = "labels_config.qs"
)