
# %% library imports
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# %% Load test obj
test = sc.read_h5ad('/Volumes/robin_work/cosmx_gray/data/samples/sample1/h5ad/sample1.h5ad')

meta = test.obs

# Define trimming thresholds
low = meta['nCount_RNA'].quantile(0.01)
high = meta['nCount_RNA'].quantile(0.99)

# Density plot
sns.kdeplot(meta['nCount_RNA'], fill=True, color="steelblue")
plt.xlabel("log1p(values)")
plt.ylabel("Density")
plt.title("Density Plot of log1p(values)")
# Add vertical lines
plt.axvline(low, color="red", linestyle="--", label=f"Low cutoff ({low:.2f})")
plt.axvline(high, color="green", linestyle="--", label=f"High cutoff ({high:.2f})")

plt.legend()
plt.show()

# %% Scaling of counts with cell area

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm

# Log transform
meta["log_area"] = np.log1p(meta["Area.um2"])
meta["log_counts"] = np.log1p(meta["nCount_RNA"])

# Fit linear model
X = sm.add_constant(meta["log_area"])
model = sm.OLS(meta["log_counts"], X).fit()

meta["fitted"] = model.fittedvalues
meta["residuals"] = meta["log_counts"] - meta["fitted"]

# Scatter with regression line
sns.scatterplot(data=meta, x="log_area", y="log_counts", alpha=0.3)
sns.lineplot(x=meta["log_area"], y=meta["fitted"], color="red", linewidth=2)
plt.xlabel("log(Cell area)")
plt.ylabel("log(Counts)")
plt.title("Scaling of counts with cell area")
plt.show()

# Threshold: top 1% residuals
threshold = meta["residuals"].quantile(0.99)
meta["flag_high"] = meta["residuals"] > threshold

# Scatter with flagged outliers
sns.scatterplot(
    data=meta, x="log_area", y="log_counts",
    hue="flag_high", palette={True: "red", False: "grey"},
    alpha=0.5, s=1
)

# Add regression line
sns.lineplot(x=meta["log_area"], y=meta["fitted"], color="blue", linewidth=2)

plt.xlabel("log(Cell area)")
plt.ylabel("log(Counts)")
plt.title("Cells with excess counts relative to size")
plt.legend(title="Flagged outlier")
plt.show()


# %%

import seaborn as sns
import matplotlib.pyplot as plt

giant_cell = meta['Area.um2'].quantile(0.1)
meta['giant_cell'] = meta['Area.um2'] < giant_cell

# Replace 'X' and 'Y' with your coordinate column names
sns.scatterplot(
    data=meta,
    x="CenterX_global_px", y="CenterY_global_px",
    hue="giant_cell",
    palette={True: "red", False: "grey"},
    alpha=0.6, s=1, markers='.'
)

plt.gca().invert_yaxis()  # optional, depending on your coordinate system
plt.title("Spatial distribution of flagged high-count cells")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.legend(title="Flagged")
plt.show()


# %% border cells

# %% test final function
test = sc.read_h5ad('/Volumes/robin_work/cosmx_gray/data/samples/sample1/h5ad/sample1.h5ad')

run_qc(test)

# Check spatially and modify rel_width accordingly
sc.pl.embedding(
    test[test.obs['fov_name'] == '1_40'], 
    basis='spatial_fov', 
    color='qc_borders', 
    size=20,
)

test_filtered = qc_filter(test)