#!/usr/bin/env python
# coding=utf-8
import sys
from Bio import SeqIO

Genome = sys.argv[1]
minLength = int(sys.argv[2])
Output = sys.argv[3]

N = 0
big_scaffolds = []

print("Reading genome...")
for scaffold in SeqIO.parse(open(Genome), "fasta"):
    N += 1
    if len(scaffold.seq) >= minLength:
        big_scaffolds.append(scaffold)

print("Found", len(big_scaffolds), "scaffolds larger than", minLength, "bp within a total of", N)
print("Writing output...")
SeqIO.write(big_scaffolds, Output, "fasta")
