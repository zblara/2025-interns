###   ＰＤＰ３ Ｍｅｔａｇｅｎｏｍｉｃｓ Ｓｔａｔｉｓｔｉｃａｌ Ａｎａｌｙｓｉｓ Ｐｉｐｅｌｉｎｅ ###
#
# This python script is optimized for the analysis of diversity statistics of filtered
# reads produced by the pdp3_mgs pipeline.
#
# The following are required as format for input data:
#   - .csv ot .tsv file
#   - Data has to be organized into a matrix containing all samples as columns and
#     the relative abundance of organisms (normalized relative counts) as rows.
#   - Tools such as Bracken (https://github.com/jenniferlu717/Bracken) + Pavian
#     (https://github.com/fbreitwieser/pavian) can be used for parsing data from Kraken 2
#
#     ©Z.B. Lara, 2024

pip install --quiet  scikit-bio seaborn matplotlib scikit-posthocs

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np  # For log transformation if needed
import pandas as pd

from google.colab import drive
drive.mount('/content/drive')

# Read the input data
data = pd.read_csv("/content/drive/My Drive/mgs_all_normalized_transposed.csv", index_col="Sample")
metadata = pd.read_csv ("/content/drive/My Drive/mgs_all_metadata.csv", index_col="Sample")
print(data)
print(metadata)

# Log transform the data
data_log = np.log1p(data.T)  # Use log1p to handle potential zeros

print(data_log)

# Calculate average abundance of each organism
average_abundance = data_log.mean(axis=1)

# Sort organisms by average abundance in descending order
sorted_organisms = average_abundance.sort_values(ascending=False)

# Select the top (y) most abundant organisms
top_organisms = sorted_organisms.head(50)

# Create a heatmap
from matplotlib.colors import LogNorm # Import LogNorm
plt.figure(figsize=(8, 10))
sns.heatmap(data_log.loc[top_organisms.index], cmap="YlGnBu", norm=LogNorm(), yticklabels=top_organisms.index)
plt.title("Most Abundant Genera Found in Sequenced Mouse Fecal Samples")
plt.xlabel("Sample")
plt.ylabel("Genus")
plt.show()

from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity

# Calculate Bray-Curtis distance matrix
distance_matrix = beta_diversity("braycurtis", data_log.T, ids=data_log.columns)
print(distance_matrix)

# Perform PCoA
pcoa_results = pcoa(distance_matrix)

# Visualize the PCoA results (example)
fig = pcoa_results.plot()
plt.show()

from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('mgs_all_w4_bray_curtis_pcoa.pdf') as pdf:
# Page 1: Bray-Curtis Distance Matrix
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(distance_matrix.data, cmap='viridis')
    ax.set_title('Bray-Curtis Distance Matrix')
    fig.colorbar(im)

    # Add values to cells
    for i in range(distance_matrix.shape[0]):
        for j in range(distance_matrix.shape[1]):
            text = ax.text(j, i, round(distance_matrix.data[i, j], 2),
                           ha="center", va="center", color="w")

    pdf.savefig(fig)
    plt.close(fig)

    # Page 2: PCoA Plot
    fig = pcoa_results.plot()  # Use the existing PCoA plot
    ax = fig.gca()  # Get the current axes
    ax.set_title('PCoA Plot')
    pdf.savefig(fig)  # Save the figure to the PDF
    plt.close(fig)  # Close the figure

# Extract data from the distance matrix
distance_array = distance_matrix.data
sample_ids = distance_matrix.ids

# Create a DataFrame from the results
distance_df = pd.DataFrame(distance_matrix, index=sample_ids, columns=sample_ids)

# Save results to a CSV file
distance_df.to_csv('mgs_all_w4_distance_matrix.csv')

# Get the coordinates of the samples in the PCoA space
sample_coordinates = pcoa_results.samples

# Print the coordinates
print(sample_coordinates)

# Get the proportion of variance explained by each principal coordinate
proportion_explained = pcoa_results.proportion_explained

# Print the proportions
print(proportion_explained)

# Extract coordinates and group information
pcoa_df = pcoa_results.samples
pcoa_df['Group'] = metadata['Group']  # Add group information from your metadata

# Group the data by 'Group'
grouped_pcoa = pcoa_df.groupby('Group')

# Calculate mean coordinates for each group
group_means = grouped_pcoa.mean()

# Calculate standard deviation of coordinates for each group
group_stds = grouped_pcoa.std()

# Save the grouped PCoA data to a CSV file
pcoa_df.to_csv('pcoa_results_by_group.csv')

# If you calculated statistics, save those as well
group_means.to_csv('pcoa_group_means.csv')
group_stds.to_csv('pcoa_group_stds.csv')

from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('mgs_w4_group_pcoa_results.pdf') as pdf:
    for group in metadata['Group'].unique(): # Iterate over unique groups in metadata
        # Subset data and metadata for the current group
        group_data = data_log.loc[:, metadata['Group'] == group]

        # Calculate Bray-Curtis distance matrix for the group
        group_distance_matrix = beta_diversity("braycurtis", group_data.T, ids=group_data.columns)

        # Perform PCoA for the group
        group_pcoa_results = pcoa(group_distance_matrix)

        # Plot PCoA for the group
        fig = group_pcoa_results.plot()
        plt.title(f"PCoA Plot for Group {group}")
        plt.show()
        print("\n")

        # Save the figure to the PDF
        pdf.savefig(fig)
        plt.close(fig) # Close the figure to avoid memory issues

        # Print proportion of variance explained for the group (also save to PDF)
        plt.figure() # Create a new figure for the text
        plt.axis('off')  # Hide axes
        plt.text(0.5, 0.5, f"Proportion of variance explained for group {group}:\n\n{group_pcoa_results.proportion_explained}",
                  ha='center', va='center', fontsize=12)
        pdf.savefig()
        plt.close()

from sklearn.decomposition import PCA

pca = PCA(n_components=3)  # Choose the number of components you want to keep
pca.fit(data_log.T)

principal_components = pca.transform(data)

plt.scatter(principal_components[:, 0], principal_components[:, 1])
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA Plot')
plt.show()

# Get unique groups and assign colors
groups = metadata['Group'].unique()
colors = plt.cm.get_cmap('tab10', len(groups))  # Choose a colormap with enough colors

# Create the PCoA plot
fig, ax = plt.subplots()

# Plot each group separately
for i, group in enumerate(groups):
    group_indices = metadata[metadata['Group'] == group].index
    ax.scatter(pcoa_results.samples.loc[group_indices, 'PC1'],
               pcoa_results.samples.loc[group_indices, 'PC2'],
               color=colors(i), label=group)

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('PCoA Plot Colored by Group')
ax.legend()
plt.show()

from matplotlib.backends.backend_pdf import PdfPages

# Get unique groups
groups = metadata['Group'].unique()

# Determine overall x and y limits across all groups
x_min, x_max = float('inf'), float('-inf')
y_min, y_max = float('inf'), float('-inf')
for group in groups:
    group_indices = metadata[metadata['Group'] == group].index
    x_vals = pcoa_results.samples.loc[group_indices, 'PC1']
    y_vals = pcoa_results.samples.loc[group_indices, 'PC2']
    x_min = min(x_min, x_vals.min())
    x_max = max(x_max, x_vals.max())
    y_min = min(y_min, y_vals.min())
    y_max = max(y_max, y_vals.max())

# Create a plot for each group with consistent scaling
with PdfPages('mgs_all_w4_pcoa_plots_by_group.pdf') as pdf:
    for group in groups:
        group_indices = metadata[metadata['Group'] == group].index
        fig, ax = plt.subplots()
        ax.scatter(pcoa_results.samples.loc[group_indices, 'PC1'],
                   pcoa_results.samples.loc[group_indices, 'PC2'])
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(f'PCoA Plot for Group {group}')
        ax.set_xlim(x_min, x_max)  # Set consistent x-axis limits
        ax.set_ylim(y_min, y_max)  # Set consistent y-axis limits

        pdf.savefig(fig)  # Save the current figure to the PDF
        plt.close(fig)  # Close the figure to avoid memory issues

from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import numpy as np

# Assuming 'distance_matrix' is your pre-calculated distance matrix and
# 'metadata' is your metadata DataFrame with a 'Group' column
mds = MDS(n_components=2, dissimilarity='precomputed', metric=False)

# Convert DistanceMatrix to a NumPy array
distance_array = distance_matrix.data
nmds_results = mds.fit_transform(distance_array)

# Get unique groups and assign colors
groups = metadata['Group'].unique()
colors = plt.cm.get_cmap('tab10', len(groups))  # Choose a colormap with enough colors

# Visualize the NMDS results with colors
for i, group in enumerate(groups):
    # Get boolean mask for the current group
    group_mask = metadata['Group'] == group
    plt.scatter(nmds_results[group_mask, 0],
                nmds_results[group_mask, 1],
                color=colors(i), label=group)

plt.xlabel('NMDS1')
plt.ylabel('NMDS2')
plt.title('NMDS Plot')
plt.legend()  # Add a legend to identify groups
plt.show()

from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Assuming 'distance_matrix' is your pre-calculated distance matrix and
# 'metadata' is your metadata DataFrame with a 'Group' column
mds = MDS(n_components=2, dissimilarity='precomputed', metric=False)

# Convert DistanceMatrix to a NumPy array
distance_array = distance_matrix.data
nmds_results = mds.fit_transform(distance_array)

# Get unique groups and assign colors
groups = metadata['Group'].unique()
colors = plt.cm.get_cmap('tab10', len(groups))  # Choose a colormap with enough colors

# Determine overall x and y limits across all groups
x_min, x_max = float('inf'), float('-inf')
y_min, y_max = float('inf'), float('-inf')
for group in groups:
    group_mask = metadata['Group'] == group
    x_vals = nmds_results[group_mask, 0]
    y_vals = nmds_results[group_mask, 1]
    x_min = min(x_min, x_vals.min())
    x_max = max(x_max, x_vals.max())
    y_min = min(y_min, y_vals.min())
    y_max = max(y_max, y_vals.max())

# Create a PDF file to save the plots
with PdfPages('mgs_all_w4_nmds_plots_by_group.pdf') as pdf:
    # Create separate plots for each group
    for i, group in enumerate(groups):
        # Get boolean mask for the current group
        group_mask = metadata['Group'] == group

        # Create a new figure for each group
        fig, ax = plt.subplots()

        ax.scatter(nmds_results[group_mask, 0],
                   nmds_results[group_mask, 1],
                   color=colors(i), label=group)

        ax.set_xlabel('NMDS1')
        ax.set_ylabel('NMDS2')
        ax.set_xlim(x_min, x_max)  # Set consistent x-axis limits
        ax.set_ylim(y_min, y_max)  # Set consistent y-axis limits
        ax.set_title(f'NMDS Plot for Group {group}')
        ax.legend()  # Add a legend to identify groups

        # Save the current figure to the PDF
        pdf.savefig(fig)

        # Close the figure to avoid memory issues
        plt.close(fig)

from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Create a PdfPages object
with PdfPages('mgs_all_w4_ordination_plots.pdf') as pdf:
  # Save the PCA plot to the PDF
  plt.figure(figsize=(8, 6))
  plt.scatter(principal_components[:, 0], principal_components[:, 1])
  plt.xlabel('Principal Component 1')
  plt.ylabel('Principal Component 2')
  plt.title('PCA Plot')
  pdf.savefig()
  plt.close()

  # Save the PCoA plot to the PDF
  fig, ax = plt.subplots()
  for i, group in enumerate(groups):
    group_indices = metadata[metadata['Group'] == group].index
    ax.scatter(pcoa_results.samples.loc[group_indices, 'PC1'],
               pcoa_results.samples.loc[group_indices, 'PC2'],
               color=colors(i), label=group)
  ax.set_xlabel('PC1')
  ax.set_ylabel('PC2')
  ax.set_title('PCoA Plot Colored by Group')
  ax.legend()
  pdf.savefig()
  plt.close()

  # Save the NMDS plot to the PDF
  # Assuming 'distance_matrix' is your pre-calculated distance matrix and
  # 'metadata' is your metadata DataFrame with a 'Group' column
  mds = MDS(n_components=2, dissimilarity='precomputed', metric=False)

  # Convert DistanceMatrix to a NumPy array
  distance_array = distance_matrix.data
  nmds_results = mds.fit_transform(distance_array)

  # Get unique groups and assign colors
  groups = metadata['Group'].unique()
  colors = plt.cm.get_cmap('tab10', len(groups))  # Choose a colormap with enough colors

  # Visualize the NMDS results with colors
  plt.figure() # Create a new figure for the NMDS plot
  for i, group in enumerate(groups):
      # Get boolean mask for the current group
      group_mask = metadata['Group'] == group
      plt.scatter(nmds_results[group_mask, 0],
                  nmds_results[group_mask, 1],
                  color=colors(i), label=group)

  plt.xlabel('NMDS1')
  plt.ylabel('NMDS2')
  plt.title('NMDS Plot')
  plt.legend()  # Add a legend to identify groups
  pdf.savefig() # Save the figure
  plt.close() # Close the figure

# Perform PERMANOVA
from skbio.stats.distance import DistanceMatrix, permanova
from statsmodels.sandbox.stats.multicomp import multipletests
import pandas as pd

ndata_transposed = data_log.T
permanova_results = permanova(distance_matrix, metadata, column='Group')

print(permanova_results)

groups = metadata['Group'].unique()

p_values = []
comparisons = []

for i in range(len(groups)):
    for j in range(i + 1, len(groups)):
        group1 = groups[i]
        group2 = groups[j]
        subset_metadata = metadata[(metadata['Group'] == group1) | (metadata['Group'] == group2)]
        subset_distance_matrix = distance_matrix.filter(subset_metadata.index)

        result = permanova(subset_distance_matrix, subset_metadata, column='Group')
        p_values.append(result['p-value'])  # Store the p-value
        comparisons.append(f"{group1} vs {group2}")  # Track comparisons for output

# Perform Bonferroni correction
reject, corrected_p_values, _, _ = multipletests(p_values, method='bonferroni')

# Print corrected p-values
for i, p_val in enumerate(corrected_p_values):
    print(f"Corrected p-value for comparison {comparisons[i]}: {p_val}")

from skbio.diversity import alpha_diversity

# Choose your metrics
metrics = ['observed_otus', 'shannon', 'simpson', 'chao1', 'fisher_alpha']

# Calculate alpha diversity for each metric
alpha_diversity_results = {}
for metric in metrics:
    alpha_diversity_results[metric] = alpha_diversity(metric, data.values, ids=data.T.columns)

# Create a pandas DataFrame from the results
alpha_diversity_df = pd.DataFrame(alpha_diversity_results)
print(alpha_diversity_df)

alpha_diversity_df.to_csv('alpha_diversity_results.csv')

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Create a PDF file to save the plots
with PdfPages('mgs_all_w4_alpha_diversity_barplots.pdf') as pdf:
    # Iterate over each alpha diversity metric
    for metric in alpha_diversity_df.columns:
        # Create a figure and axes for the plot
        fig, ax = plt.subplots()

        # Create the bar plot with standard error
        sns.barplot(x='Group', y=metric, data=metadata.join(alpha_diversity_df), ax=ax, errorbar='se')

        # Set the title of the plot
        ax.set_title(f'Alpha Diversity ({metric}) by Group')

        # Rotate x-axis labels if needed
        plt.xticks(rotation=45, ha='right')

        # Save the figure to the PDF
        pdf.savefig(fig)
        plt.close(fig)  # Close the figure to avoid memory issues

# Combine results into a single DataFrame
alpha_div_df = pd.DataFrame(alpha_diversity_results)
alpha_div_df.index.name = 'Sample'

# Scale each metric to a common range (0-1)
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
alpha_div_df[metrics] = scaler.fit_transform(alpha_div_df[metrics])

# Merge with metadata
alpha_div_df = pd.merge(alpha_div_df, metadata, left_index=True, right_index=True)

# Melt the DataFrame to long format for plotting
alpha_div_melt = pd.melt(alpha_div_df, id_vars=['Group'], value_vars=metrics,
                        var_name='Metric', value_name='Diversity')

# Plot grouped column plot with dot patterns and shading
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.figure(figsize=(10, 6))
markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd']

# Use a pastel color palette
palette = sns.color_palette("pastel", n_colors=len(groups))

ax = sns.barplot(x='Metric', y='Diversity', hue='Group', data=alpha_div_melt,
                 errorbar='se', palette=palette, edgecolor='black')

# Create a PDF file to save the plots
with PdfPages('mgs_all_w4_alpha_diversity_metrics_plot.pdf') as pdf:
    # Add dot patterns to bars
    for i, bar in enumerate(ax.patches):
        group_idx = i // len(metrics)
        # Ensure group_idx is within the valid range of palette indices
        group_idx = group_idx % len(palette)
        bar.set_facecolor(palette[group_idx])  # Set bar color
        bar.set_edgecolor('black')  # Add black outline for better visibility

        # Calculate dot positions within the bar
        bar_width = bar.get_width()
        bar_height = bar.get_height()
        x_center = bar.get_x() + bar_width / 2
        y_center = bar.get_y() + bar_height / 2

        # Scatter plot for dots
        ax.scatter(x_center, y_center, marker=markers[group_idx], color='white', s=50, zorder=3)  # Adjust dot size as needed

    # Create custom legend elements with markers
    legend_elements = []
    for i, group in enumerate(groups):
        legend_elements.append(
            plt.Line2D([0], [0], marker=markers[i], color='w', label=group,
                      markerfacecolor=palette[i], markersize=10, ls='None'))  # Adjust markersize as needed

    plt.title('Alpha Diversity Metrics by Group')
    plt.xlabel('Metric')
    plt.ylabel('Diversity (Scaled 0-1)')
    plt.xticks(rotation=45, ha='right')
    plt.legend(handles=legend_elements, title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    pdf.savefig()
    plt.show()
    plt.close(fig)  # Close the figure to avoid memory issues

from scipy import stats

# Iterate over each alpha diversity metric
for metric in alpha_diversity_df.columns:
    # Group data by the 'Group' column in metadata
    groups = [alpha_diversity_df.loc[metadata[metadata['Group'] == g].index, metric] for g in metadata['Group'].unique()]

    # Perform Kruskal-Wallis test
    statistic, p_value = stats.kruskal(*groups)

    print(f"Kruskal-Wallis test for {metric}:")
    print(f"H statistic: {statistic}")
    print(f"P-value: {p_value}")

# Kruskal-Wallis between groups
import scikit_posthocs as sp
from statsmodels.sandbox.stats.multicomp import multipletests

for metric in alpha_diversity_df.columns:
    # Group data by the 'Group' column in metadata
    groups = [alpha_diversity_df.loc[metadata[metadata['Group'] == g].index, metric].values for g in metadata['Group'].unique()]

    # Perform Dunn's test
    p_values = sp.posthoc_dunn(groups, p_adjust='bonferroni') # Use Bonferroni correction for initial pairwise comparisons

    # Apply Benjamini-Hochberg correction to the p-values from Dunn's test
    # Extract the NumPy array from the DataFrame before flattening
    reject, corrected_p_values, _, _ = multipletests(p_values.values.ravel(), method='fdr_bh') # 'fdr_bh' is Benjamini-Hochberg

    # Reshape corrected p-values back to a matrix
    corrected_p_values = corrected_p_values.reshape(p_values.shape)

    print(f"Dunn's Test with Benjamini-Hochberg correction for {metric}:")
    print(corrected_p_values)

#Kruskal Wallis Organisms
import pandas as pd
import numpy as np
import scipy.stats as stats
import scikit_posthocs as sp
from statsmodels.sandbox.stats.multicomp import multipletests

# Choose the organism
organism_name = "Lactiplantibacillus"  # Replace with the actual organism name

# Extract data for the organism
data_extr = data.T
organism_data = data_extr.loc[organism_name]

# Group data by 'Group' in metadata
groups = [organism_data[metadata['Group'] == group].values for group in metadata['Group'].unique()]
groups_mean = [np.mean(organism_data[metadata['Group'] == group]) for group in metadata['Group'].unique()]
print(f"{organism_name}:\n {groups_mean}")

# Perform Kruskal-Wallis test
statistic, p_value = stats.kruskal(*groups)
print(f"Kruskal-Wallis statistic: {statistic}")
print(f"P-value: {p_value}")

# Perform Dunn's test with Bonferroni correction
dunn_results = sp.posthoc_dunn(groups, p_adjust='bonferroni')
print("Dunn's Test with Bonferroni correction:")
print(dunn_results)

# Flatten the p-value matrix from Dunn's test
flat_p_values = dunn_results.values.ravel()

# Apply Benjamini-Hochberg correction to the flattened p-values
reject, corrected_p_values, _, _ = multipletests(flat_p_values, method='fdr_bh')

# Reshape corrected p-values back to the original matrix shape
corrected_p_values = corrected_p_values.reshape(dunn_results.shape)

# Create a DataFrame for the corrected p-values for better readability
corrected_p_values_df = pd.DataFrame(corrected_p_values, index=dunn_results.index, columns=dunn_results.columns)

print("Dunn's Test with Benjamini-Hochberg correction:")
print(corrected_p_values_df)

import scipy.stats as stats
import pandas as pd
import numpy as np

# Initialize a dictionary to store results
results = {}
data_extr = data.T

# Iterate through each organism
for organism_name in data_extr.index:
    # Extract data for the organism
    organism_data = data_extr.loc[organism_name]

    # Group data by 'Group' in metadata
    groups = [organism_data[metadata['Group'] == group].values for group in metadata['Group'].unique()]
    groups_mean = [np.mean(organism_data[metadata['Group'] == group]) for group in metadata['Group'].unique()]
    print(groups_mean)

    # Check for zero variance before performing the test
    if np.any([np.var(g) != 0 for g in groups]):  # Check if any group has non-zero variance
        # Perform Kruskal-Wallis test
        statistic, p_value = stats.kruskal(*groups)
        print(f"{organism_name}: Kruskal-Wallis statistic: {statistic}")
        print(f"P-value: {p_value}")

        # Store results
        results[organism_name] = {'statistic': statistic, 'p_value': p_value}
    else:
        print(f"Skipping {organism_name} due to zero variance within groups.")

# Create a DataFrame from the results
results_df = pd.DataFrame.from_dict(results, orient='index')

# Save results to a CSV file
results_df.to_csv('mgs_organisms_w4_kruskal_wallis_results.csv')

import pandas as pd
import numpy as np
from scipy import stats
import scikit_posthocs as sp
from statsmodels.sandbox.stats.multicomp import multipletests

# Initialize a list to store significant organisms and their Dunn's test results
significant_organisms = []
results = {}
data_extr = data.T

# Iterate through each organism
for organism_name in data_extr.index:
    # Extract data for the organism
    organism_data = data_extr.loc[organism_name]

    # Group data by 'Group' in metadata
    groups = [organism_data[metadata['Group'] == group].values for group in metadata['Group'].unique()]
    groups_mean = [np.mean(organism_data[metadata['Group'] == group]) for group in metadata['Group'].unique()]
    print(groups_mean)

    # Check for zero variance before performing the test
    if np.any([np.var(g) != 0 for g in groups]):  # Check if any group has non-zero variance
        # Perform Kruskal-Wallis test on the original data, not the means
        statistic, p_value = stats.kruskal(*groups)  # Use 'groups' here
        print(f"{organism_name}: Kruskal-Wallis statistic: {statistic}")
        print(f"P-value: {p_value}")

        # Store the organism if p-value is significant
        if p_value < 0.05:
            significant_organisms.append(organism_name)

            # Perform Dunn's test with Bonferroni correction
            dunn_bonferroni = sp.posthoc_dunn(groups, p_adjust='bonferroni')

            # Perform Dunn's test and apply Benjamini-Hochberg correction
            dunn_p_values = sp.posthoc_dunn(groups)
            reject, corrected_p_values, _, _ = multipletests(dunn_p_values.values.ravel(), method='fdr_bh')
            dunn_bh = pd.DataFrame(corrected_p_values.reshape(dunn_p_values.shape),
                                  index=dunn_p_values.index, columns=dunn_p_values.columns)

            # Store results
            results[organism_name] = {'Kruskal-Wallis p-value': p_value,
                                       'Dunn Bonferroni': dunn_bonferroni,
                                       'Dunn Benjamini-Hochberg': dunn_bh}
    else:
        print(f"Skipping {organism_name} due to zero variance within groups.")

# Save significant organisms and Dunn's test results to a text file
with open('mgs_all_w4_significant_organisms_with_dunns.txt', 'w') as f:
    for organism_name in significant_organisms:
        f.write(f"Organism: {organism_name}\n")
        f.write(f"Kruskal-Wallis p-value: {results[organism_name]['Kruskal-Wallis p-value']}\n")

        f.write("\nDunn's Test with Bonferroni correction:\n")
        f.write(results[organism_name]['Dunn Bonferroni'].to_string())

        f.write("\nDunn's Test with Benjamini-Hochberg correction:\n")
        f.write(results[organism_name]['Dunn Benjamini-Hochberg'].to_string())

        f.write("\n--------------------\n")

import pandas as pd

# Load the dataframe.
df = pd.read_csv("/content/drive/My Drive/mgs_all_normalized_transposed.csv", index_col='Sample')

# Extract group names from sample IDs.
df['GroupID'] = df.index.str[0]

# Calculate the average for each group.
df = df.groupby('GroupID').mean()

# Save the averaged data to a new CSV file.
df.to_csv('averaged_data.csv')
avdata = df.T
avdata.to_csv('mgs_averaged_data.csv')

import pandas as pd
from scipy.stats import spearmanr

# Verify the file name and path
abundance_data = pd.read_csv("/content/drive/My Drive/averaged_data.csv", index_col=0)
serum_metadata = pd.read_csv("/content/drive/My Drive/mgs_serum_metadata.csv", index_col="Sample")

results = {}
for organism in abundance_data.columns:
    results[organism] = {}
    for lipid in ['CHOLE', 'TG', 'HDL', 'LDL']:
        correlation, p_value = spearmanr(abundance_data[organism], serum_metadata[lipid])
        results[organism][lipid] = {'correlation': correlation, 'p_value': p_value}

results_df = pd.DataFrame.from_dict({(i,j): results[i][j]
                               for i in results.keys()
                               for j in results[i].keys()},
                            orient='index')
results_df.index = pd.MultiIndex.from_tuples(results_df.index)

# Save the results to a txt file
results_df.to_csv('mgs_serum_correlation_results.csv')