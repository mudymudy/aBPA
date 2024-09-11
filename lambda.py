#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 13:58:15 2024

@author: bruno
"""

import sys
import pandas as pd
import os



data = sys.argv[1]
base_filename = os.path.basename(data)
output_file = f"{base_filename}_final.csv"

with open(data, 'r') as file:
        summary = pd.read_csv(data, sep='\t')


updated_series = []

for index, row in summary.iterrows():
    if row['completeness'] >= 50 and row['normalizedCoverage'] >= 0.5:
        updated_series.append(1)
    else:
        updated_series.append(0)


final_df = pd.DataFrame({'Gene': summary['Gene'], 'presenceAbsence': updated_series})
final_df.to_csv(output_file, index=False)
