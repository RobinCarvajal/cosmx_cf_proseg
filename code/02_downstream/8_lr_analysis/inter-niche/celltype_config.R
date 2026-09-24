library(whereami)

wd <- dirname(thisfile())
setwd(wd)

## Map spatial niches to broader biological compartments
map.celltype_group <- c(
  # STROMAL
  Stromal_Interstitium   = "STROMA",
  Alveolar_Interstitium  = "STROMA",

  # LYMPHOID
  Plasma_Rich            = "LYM",
  BALT                   = "LYM",

  # MYELOID
  Neutrophil_Rich        = "MYE",

  # AIRWAY EPITHELIUM
  Basal_Airway           = "AIR_EPI",
  Luminal_Airway         = "AIR_EPI",

  # ALVEOLAR EPITHELIUM
  Alveolar_Epithelium    = "ALV_EPI"
)

## Desired group ordering
levels.group <- c(
  "ALV_EPI",
  "AIR_EPI",
  "MYE",
  "LYM",
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