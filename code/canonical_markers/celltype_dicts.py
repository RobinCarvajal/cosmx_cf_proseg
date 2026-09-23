epithelial_markers = {
    'basal': ['KRT5', 'KRT14', 'KRT15', 'KRT17', 'IL33'],
    'secretory': ['SCGB1A1', 'SCGB3A1'],
    'ciliated': ['FOXJ1', 'TUBB4B', 'SAA1'],
    'trSCs' : ['SCGB1A1', 'SCGB3A1', 'FOXJ1', 'TUBB4B', 'SAA1'],
    'AT1': ['RTKN2', 'COL4A3', 'FSTL3', 'AGER'], 
    'AT2': ['NKX2-1', 'ETV5', 'NAPSA', 'LAMP3', 'CD36', 'LPCAT1', 'SFTPD'],
    'hillock': ['KRT13', 'KRT6A', 'DSG3', 'SERPINB2'],
}

lymphoid_markers = {
    'T': ['CD3D','CD3G','CD3E','TRBC2','CD8A','CD8B'],
    'B': ['MS4A1','CD19','CD74','CD79A','IGHD','IGHM'],
    'B_GERM': ['MZB1','IGLL1','MKI67'],
    'Plasma': ['IGKC','JCHAIN','IGHG1/2','IGHA1'],
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


stroma_markers = {
    'fibroblast': ['PDGFRB','LUM','COL6A2','COL1A1','COL1A2'],
    'smooth_muscle': ['ACTA2','MYL9','RGS5','NOX4','TAGLN','CALD1']
}

endothelial_markers = {
    'EC': ['PECAM1','CD93','VWF','EMCN','FLT1','ID3','MCAM','CLDN5']
}

all_markers = epithelial_markers | \
    lymphoid_markers |  \
    myeloid_markers | \
    stroma_markers | \
    endothelial_markers


# make all_markers a JSON file 
import json
# get the parent directory of this file
import os
parent_dir = '/data/cosmx_gray/code/canonical_markers'
with open(os.path.join(parent_dir, 'all_markers.json'), 'w') as f:
    json.dump(all_markers, f, indent=4)

