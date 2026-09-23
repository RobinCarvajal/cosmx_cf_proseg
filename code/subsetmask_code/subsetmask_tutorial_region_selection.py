# -*- coding: utf-8 -*-
"""
Created on Sat Jun 21 11:11:13 2025

@author: robin
"""

# %% Run this before to enable qt (needed for napari)
%gui qt

# %% (SAMPLE LEVEL) create the Image object and store an image 


import pandas as pd
import subsetmask as sm

# Set the working directory
import os
os.chdir("/media/robin/robin_work/cosmx_gray")

# load samples_meta with samples 
samples_meta = pd.read_csv("data/metadata_files/slide_1_metadata_with_samples.csv")

# test
sample_1 = sm.Image(metadata_df=samples_meta, 
                  images_col='samples',
                  image_name='sample_1',
                  x_col="x_slide_mm",
                  y_col="y_slide_mm")

# addd store image method
sample_1.store_image(group_col=None, figsize=(30, 25), dpi=300, dot_size=5, cmap='tab20')

## %%
sample_1.draw_labels()

#%% save sample_1 annotations
# without closing napari we run this line to save the annotations
sample_1.save_annotations(ann_col_name="regions")

# %% we change our annotations for something more meaningful
map = {
    "label_1": "region_1",
    "label_2": "region_2",
    "label_3": "region_3"
}

sample_1.rename_annotations(ann_col_name="regions", mapping_dict=map)

# %% we extract the annotated df from the object
regions_meta = sample_1.annotated_metadata_df

regions_meta


# %% Checking that it worked

# scater plot using matplot lib
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(7, 7))
sns.set_style("white")
sns.scatterplot(
    data=regions_meta,
    x="x_slide_mm",
    y="y_slide_mm",
    hue="regions",        # or "library_id", etc.
    s=3,
    alpha=0.7,
    palette="tab20"      # or any seaborn/matplotlib palette
)
plt.axis("equal")
plt.title("Centroids colored by region")
plt.show()

# %% Save the annotations metadata

