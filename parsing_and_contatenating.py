#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno
"""

from Bio import SeqIO
import os
import sys


path = sys.argv[1]


# Path to the directory containing downloaded GenBank files
genbank_dir = "path"
output_fasta_file = "clustered_sequences.fasta"

all_genes = {}

# Function to process GenBank file and extract gene information
def process_genbank(gb_file):
    with open(gb_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    gene_name = feature.qualifiers.get("gene", ["unknown_gene"])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown_locus"])[0]
                    gene_sequence = feature.location.extract(record).seq
                    all_genes[locus_tag] = {
                        "gene_name": gene_name,
                        "sequence": str(gene_sequence)
                    }

# Iterate over all GenBank files in the directory
for gb_file in os.listdir(genbank_dir):
    if gb_file.endswith(".gbff"):
        process_genbank(os.path.join(genbank_dir, gb_file))

# Write the collected gene information to a single FASTA file
with open(output_fasta_file, "w") as output_handle:
    for locus_tag, gene_info in all_genes.items():
        output_handle.write(f">{gene_info['gene_name']}|{locus_tag}\n{gene_info['sequence']}\n")

print(f"All gene sequences have been written to {output_fasta_file}")
