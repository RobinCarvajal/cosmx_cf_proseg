epithelial_markers = {
    'basal': ['KRT5', 'KRT14', 'KRT6C', 'KRT15'],
    'parabasal': ['KRT4', 'KRT13', 'KRT16', 'KRT23'],
    'club': ['SCGB1A1', 'SCGB3A1'],
    'goblet': ['SPDEF'],
    'ciliated': ['FOXJ1', 'TUBB4B', 'SAA1', 'MYB', 'RFX2', 'RFX3', 'TP73'],
    'trSCs' : ['SCGB1A1', 'SCGB3A1', 'FOXJ1', 'TUBB4B', 'SAA1'],
    'deuterosomal': ['FOXJ1','CDC20B', 'CCNO'],
    'tuft': ['POU2F2', 'POU2F3'],
    'NE': ['NGF'],
    'hillock': ['KRT13', 'KRT6A', 'DSG3', 'SERPINB2'],
    'hillock-club': ['KRT13', 'SCGB1A1'],
    'AT1': ['RTKN2', 'COL4A3', 'FSTL3', 'AGER'], 
    'AT2': ['NKX2-1', 'ETV5', 'NAPSA', 'LAMP3', 'CD36', 'LPCAT1', 'SFTPD'],
}



# make all_markers a JSON file 
import json
# get the parent directory of this file
import os
parent_dir = '/data/cosmx_gray/code/canonical_markers'
with open(os.path.join(parent_dir, 'epithelial_markers.json'), 'w') as f:
    json.dump(epithelial_markers, f, indent=4)

