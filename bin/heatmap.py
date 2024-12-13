#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 10:56:18 2024

@author: bruno
"""
import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
sys.setrecursionlimit(20000)
import fastcluster

'''
data = the final matrix
names = unique_names. This is so we can detect the actual samples in the dataset and do things
'''
matplotlib.use('Agg')

data=sys.argv[1]
namesfile=sys.argv[2]

with open(namesfile, 'r') as f:
    names = f.read().split()  # Read all names and split by whitespace

print("Sample names loaded: ", names)

with open(data, 'r') as file:
        matrix = pd.read_csv(data, sep='\t')
        print("matrix.tab file loaded successfully into python.")
        print(matrix.head())

matrix = matrix.set_index(matrix.columns[0]) 


"""
I need to do:
                - Select list of genes that are present in at least one ancient strain.
                - Select list of genes that are present in every sample.
                - Select list of genes that excludes genes that are present in every sample.
                - Get list of sample names from the clustered heatmap, so we can know which samples are closer to ancient strains.
"""

#Select list of genes that excludes genes that are present in every sample.
mask = (matrix == 1).all(axis=1)
maskedMatrixNoUbiquitous = matrix[~mask]
maskedMatrixGenesNoUbiquitous = maskedMatrixNoUbiquitous.index
maskedMatrixGenesNoUbiquitous_series = pd.Series(maskedMatrixGenesNoUbiquitous, name="Genes")
maskedMatrixGenesNoUbiquitous_series.to_csv("maskedMatrixGenesNoUbiquitous.txt", sep="\t", index=False, header=False)



#Select list of genes that are present in at least one ancient strain.
matrixOnlyAncient = matrix[names]
maskOnlyAncient = (matrixOnlyAncient == 0).all(axis=1)
maskedOnlyAncient = matrix[~maskOnlyAncient]
maskedMatrixGenesOnlyAncient = maskedOnlyAncient.index
maskedMatrixGenesOnlyAncient_series = pd.Series(maskedMatrixGenesOnlyAncient, name="Genes")
maskedMatrixGenesOnlyAncient_series.to_csv("maskedMatrixGenesOnlyAncient.txt", sep="\t", index=False, header=False)


#Select list of genes that are present in every sample.
maskedMatrixUbiquitous =  matrix[mask]
maskedMatrixGenesUbiquitous = maskedMatrixUbiquitous.index
maskedMatrixGenesUbiquitous_series = pd.Series(maskedMatrixGenesUbiquitous, name="Genes")
maskedMatrixGenesUbiquitous_series.to_csv("maskedMatrixGenesUbiquitous.txt", sep="\t", index=False, header=False)

# 1. Get the number of samples in the matrix, which will be 100%
# 2. Iterate over each row and select genes in the matrix that are present at least 80% of the samples
# 3. Make a new dataframe with this set of genes: Plot it and export it to working directory.

num_samples = matrix.shape[1]
threshold = 0.8 * num_samples
genesAbovePercent = matrix[(matrix.sum(axis=1) >= threshold)]
genesAbovePercentIndex = genesAbovePercent.index
genesAbovePercentSeries = pd.Series(genesAbovePercentIndex, name="Genes")
genesAbovePercentSeries.to_csv("genesAbovePercentSeries.txt", sep="\t", index=False, header=False)




'''
I NEED TO IMPLEMENT THIS IN A WAY THAT IT WILL DETECT THE ACTUAL SAMPLES

# Select the last two columns
matrix_02 = matrix.iloc[:, -2:]
# Create a boolean mask where both columns are 0
mask_05 = (matrix_02.iloc[:, 0] == 0) & (matrix_02.iloc[:, 1] == 0)

# Filter the DataFrame to keep only rows where the mask is False
filtered_matrix_05 = matrix_05[~mask_05]

'''



def create_clustered_heatmap(data, names, label):
    """
    Creates a clustered heatmap with the specified formatting and customization.
    
    Parameters:
    data (pandas.DataFrame): The input data for creating the heatmap.
    
    Returns:
    None
    """
    ancient_strains = names

    plt.figure(figsize=(50, 120), dpi=1000)

    sns.set(font_scale=1)

    colors = ["orange", "darkred"]
    cmap = ListedColormap(colors)
    line_kws = {"linewidth": 4}
    
    clustered_heatmap = sns.clustermap(
        data, 
        row_cluster=True, 
        col_cluster=True, 
        cmap=cmap, 
        method="complete", 
        metric='euclidean', 
        figsize=(50, 120), 
        tree_kws=line_kws, 
        dendrogram_ratio=(0.25, 0.25)
    )

    clustered_heatmap.ax_heatmap.tick_params(
        axis='y', which='both', bottom=True, top=False, labelbottom=True, labeltop=False
    )

    clustered_heatmap.ax_heatmap.set_yticklabels([])
    
    colorbar = clustered_heatmap.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['Absence', 'Presence'])
    
    plt.setp(clustered_heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontweight='bold', fontsize=10)
    
    plt.setp(clustered_heatmap.ax_heatmap.yaxis.get_label(), fontweight='bold', fontsize=10)
    clustered_heatmap.ax_heatmap.set_ylabel("")
    
    colorbar.ax.tick_params(labelsize=50)
    
    #Get the positions of columns in the clustered heatmap
    x_ticks = clustered_heatmap.ax_heatmap.get_xticks()
    x_labels = [label.get_text() for label in clustered_heatmap.ax_heatmap.get_xticklabels()]

    #Create a mapping from column names to their positions
    column_to_position = {label: pos for label, pos in zip(x_labels, x_ticks)}
    
    def extend_positions(positions, data_shape):
        """Extend positions to include neighbors but avoid redundant lines."""
        extended = set()
        for pos in positions:
            if 0 < pos < data_shape[1]:
                extended.add(pos + 0.5)  #Add the position itself
            if pos - 1 >= 0:
                    extended.add(pos - 0.5)  #Add the left neighbor
        return sorted(extended)
    
    def plot_lines(positions, color):
        """Plot vertical lines at specific columns with adjustments."""
        unique_positions = sorted(set(positions))  #Remove duplicates
        for vline in unique_positions:
            if 0 <= vline < data.shape[1]:
                clustered_heatmap.ax_heatmap.axvline(x=vline, color=color, linewidth=3)
        print(f"Plotted lines at positions: {unique_positions}")


    #Find positions of ancient strains and extend to include neighbors
    ancient_strains_pos = extend_positions([column_to_position[col] for col in ancient_strains if col in column_to_position], data.shape)
    plot_lines(ancient_strains_pos, 'blue')
    
    ordered_columns_indices = clustered_heatmap.dendrogram_col.reordered_ind
    sampleOrder = data.columns[ordered_columns_indices]

    sampleOrder_series =  pd.Series(sampleOrder, name="Species")
    sampleOrder_series.to_csv(f"sampleOrder{label}.txt", sep="\t", index=False, header=False)
    

    #Add vertical line at position 0, NOT ANYMORE
    #clustered_heatmap.ax_heatmap.axvline(x=0, color='blue', linewidth=3)
    plt.savefig(f"presenceAbsenceMatrix{label}.png",bbox_inches='tight') 
    plt.tight_layout()
    plt.show()


create_clustered_heatmap(maskedMatrixNoUbiquitous, names, "noUbiquitous")
create_clustered_heatmap(maskedOnlyAncient, names, "onlyAncient")
create_clustered_heatmap(genesAbovePercent, names, "abovePercent")
