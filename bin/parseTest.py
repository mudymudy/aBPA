#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno
"""

from Bio import SeqIO
import sys
import os

def verify_genbank_file(gb_file):
    try:
        with open(gb_file, "r") as handle:
            # Attempt to parse the file using Biopython
            for _ in SeqIO.parse(handle, "genbank"):
                pass
        print(f"{gb_file} is a valid GenBank file.")
        return True
    except Exception as e:
        print(f"{gb_file} is not a valid GenBank file. Error: {e}")
        return False

# Directory containing GenBank files
genbank_dir = sys.argv[1]

# Iterate over all files in the directory and verify them
for gb_file in os.listdir(genbank_dir):
    if gb_file.endswith(".gbff"):
        verify_genbank_file(os.path.join(genbank_dir, gb_file))
