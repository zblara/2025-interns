#####################################################################################
#
#  ZBL Metagenomics Classification Pipeline
#
#  This bash script is optimized for the normalization of taxonomic reads.
#
#  Usage:
#	Replace input data file with your .csv file containing raw abundance counts
#   from Kraken2 and Pavian tools.
#	
#  Example:
#	Pavian reads: /media/zbl/Storage/mgs(/pavian)
#	
#  Output:
#  	1. CLR-normalized data (*_CLR_transformed_*.csv")
#
#  Date: 2025-05-02
#  Author: Z.B. Lara 
#
######################################################################################

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read input data
data = pd.read_csv("/mnt/d/Guiuan/statistics/Guiuan_taxonomy_genus.csv", index_col="name")
metadata = pd.read_csv("/mnt/d/Guiuan/statistics/Guiuan_metadata.csv", index_col="Sample")

print(data.head())
print(metadata.head())


# Transpose so samples are rows and taxa are columns
data_T = data.T

# Add a pseudo-count to avoid zeros before CLR
pseudo_count = 1e-6
data_T_pseudo = data_T + pseudo_count

# CLR Transformation
def clr(x):
    return np.log(x / np.exp(np.mean(np.log(x))))

clr_data = data_T_pseudo.apply(clr, axis=1)

# Transpose back to have taxa as rows
clr_data_T = clr_data.T

# Calculate average abundance per taxa
average_abundance = clr_data_T.mean(axis=1)
sorted_organisms = average_abundance.sort_values(ascending=False)

# Select top 50 most abundant taxa
top_organisms = sorted_organisms.head(50)
data_top = clr_data_T.loc[top_organisms.index]

# Export the full CLR-transformed matrix and top 50 taxa to CSV
clr_data_T.to_csv("/mnt/d/Guiuan/statistics/Guiuan_CLR_transformed_genus.csv")

print("CLR-normalized data exported successfully.")
