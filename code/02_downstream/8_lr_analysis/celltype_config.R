library(whereami)

wd <- dirname(thisfile())
setwd(wd)

## Map detailed cell-type labels to broader groups
map.celltype_group <- c(
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
  lymEC = "EC"

)

## Desired group ordering
levels.group <- c(
  "ALV_EPI",
  "AIR_EPI",
  "MYE",
  "LYM",
  "EC",
  "STROMA"
)

## Store both objects together
celltype_config <- list(
  map.celltype_group = map.celltype_group,
  levels.group = levels.group
)

## Save configuration
saveRDS(
  celltype_config,
  file = "celltype_config.rds"
)