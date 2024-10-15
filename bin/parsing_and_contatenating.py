#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 16:14:48 2024

@author: bruno
"""
from Bio import SeqIO
import os
import sys

genbank_dir = sys.argv[1]
output_fasta_file = "clustered_sequences.fasta"

all_genes = {}

# Function to process GenBank file and extract gene information
def process_genbank(gb_file):
    with open(gb_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                # Extract locus_tag, gene, and product information
                locus_tag_name = feature.qualifiers.get("locus_tag", ["unknown_locus_tag"])[0]
                gene_name = feature.qualifiers.get("gene", [""])[0]  # Use empty string if gene is absent
                product_name = feature.qualifiers.get("product", [""])[0]  # Use empty string if product is absent

                # Clean up the gene name and product name by replacing problematic characters
                product_name = product_name.replace("/", "_").replace("+", "_")

                # Extract gene sequence
                gene_sequence = feature.location.extract(record).seq
                all_genes[locus_tag_name] = {
                    "gene_name": gene_name,
                    "product_name": product_name,
                    "sequence": str(gene_sequence)
                }

# Iterate over all GenBank files in the directory
for gb_file in os.listdir(genbank_dir):
    if gb_file.endswith(".gbff"):
        process_genbank(os.path.join(genbank_dir, gb_file))

# Write the collected gene information to a single FASTA file
with open(output_fasta_file, "w") as output_handle:
    for locus_tag_name, gene_info in all_genes.items():
        # Write the output in the desired format
        output_handle.write(f">{locus_tag_name} ~~~{gene_info['gene_name']}~~~{gene_info['product_name']}~~~\n{gene_info['sequence']}\n")

print(f"All gene sequences have been written to {output_fasta_file}")
