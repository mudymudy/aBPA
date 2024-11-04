#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:20:12 2024

@author: bruno
"""

from Bio import AlignIO
import sys
import os



alignment = sys.argv[1]
alignmentDir = os.path.dirname(alignment)
output = os.path.join(alignmentDir, "pMauveFastaMSA.fasta")

AlignIO.convert(alignment, "mauve", output ,"fasta" )
