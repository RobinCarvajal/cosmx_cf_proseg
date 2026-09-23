
# make all_markers a JSON file 
import json
# get the parent directory of this file
import os


lymphoid_markers = {
    'T': ['CD3D','CD3G','CD3E'],
    'CD4+ T': ['CD4','IL7R'],
    'CD8+ T': ['CD8A','CD8B'],
    'NKT': ['GNLY','CD247','PTGDR','SH2D1B'],
    'B': ['MS4A1','CD19','CD74','CD79A','IGHD'],
    'PC': ['IGKC','JCHAIN','IGLC1/2','IGLL5'],
    'PC_IgG': ['IGHG1/2'],
    'PC_IgA': ['IGHA1'],
    'PC_IgM': ['IGHM'],
}

myeloid_markers = {
    'APC': ['HLA-DRA','HLA-DRB'],
    'neutrophils': ['ITGAM','CEACAM8','CSF3R','ELANE','S100A8','S100A9'],
    'mast': ['TPSAB1/2','KIT','GATA2'],
    'DC': ['ITGAM','IRF8','CD83'],
    # ITGAX(CD11C) needs to be low expression
    'macrophages': ['ITGAM','CD74','LYZ','CD68','ITGAX'],
    # CD86 and MRC1(CD206) need to be low expression 
    # added CD16 (FCGR3A/CD16a) 
    'monocytes': ['CD14','FCGR3A','CD86','CCR2','MRC1','APOE'],  
    # MARCO and MRC1(CD206) for alveolar macrophages, the others not so sure
    'alveolar_macrophages': ['MARCO','MRC1','PPARG','CD163','CD36'],
}


all_markers = lymphoid_markers |  \
    myeloid_markers 



parent_dir = '/mnt/data/project0062/proseg_data/code/canonical_markers'
with open(os.path.join(parent_dir, 'immune_markers.json'), 'w') as f:
    json.dump(all_markers, f, indent=4)

