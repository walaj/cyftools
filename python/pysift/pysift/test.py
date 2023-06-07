import scimap as sm
import anndata as ad
import pandas as pd

# load the data
data = pd.read_csv('/Users/jeremiahwala/Sorger/orion/orion_1_74/LSP10452.rad.csv')
meta = pd.read_csv('/Users/jeremiahwala/Sorger/orion/orion_1_74/LSP10452.meta.csv')

# Select columns 3, 4 (1-based), which are 2, 3 (0-based)
cols_xy = data.iloc[:, [2, 3]]

# Select columns starting from 27 (1-based), which is 26 (0-based)
cols_meta = data.iloc[:, 26:]

# Concatenate the two DataFrames along the column axis
subset_data = pd.concat([cols_xy, cols_meta], axis=1)



# Select rows starting from 27 (1-based), which is 26 (0-based)
subset_meta = meta.iloc[26:]

# Create a DataFrame with "x" and "y"
xy = pd.DataFrame(["x", "y"], columns=subset_meta.columns)

# Concatenate the two DataFrames along the row axis
subset_meta = pd.concat([xy, subset_meta], ignore_index=True)

subset_data = data.iloc[:, 26:]
adata = ad.AnnData(subset_data)

