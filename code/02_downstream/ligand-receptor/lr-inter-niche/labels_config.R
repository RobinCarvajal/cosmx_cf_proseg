library(whereami)
library(qs)

wd <- dirname(thisfile())
setwd(wd)

## Map spatial niches to broader biological compartments
label_to_group <- c(
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

# all groups
all_groups <- unique(unname(label_to_group))

## Store both objects together
labels_config <- list(
  label_to_group = label_to_group,
  all_groups = all_groups
)

## Save configuration
qsave(
  labels_config,
  file = "labels_config.qs"
)