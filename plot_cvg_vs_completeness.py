#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle

data = sys.argv[1]
completeness = sys.argv[2]
coverage = sys.argv[3]


with open(data, 'r') as file:
        summary = pd.read_csv(data, sep='\t')
        print("geneNormalizedSummary.tab file loaded successfully into python.")
        print(summary.head())

max_value = summary['normalizedGeneSimple'].max()

line_segments = [
    [(completeness, 0), (completeness, max_value)],  # Line 1: vertical line at x=50
    [(0, coverage), (100, coverage)]  # Line 2: horizontal line at y=1.5
]

plt.figure(figsize=[10,10])
ax = sns.scatterplot(data=summary, y='normalizedGeneSimple' , x='geneCompleteness' , hue='sampleID', s=7)

lines = LineCollection(line_segments, linewidths=2, colors='black', linestyle='dashed')
ax.add_collection(lines)


# Create the Rectangle patch
rec1 = Rectangle((0, 0), 50, 12, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)
rec2 = Rectangle((50, 0), 50, 0.5, linewidth=1, edgecolor='none', facecolor='grey', alpha=0.3)

# Add the rectangle to the plot
ax.add_patch(rec1)
ax.add_patch(rec2)

plt.legend(markerscale=1, fontsize=10, loc='upper left', bbox_to_anchor=(-0.01,1.1))


plt.savefig('plotCoverage_vs_Completeness.png') 
