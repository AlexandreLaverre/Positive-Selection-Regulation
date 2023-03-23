#!/usr/bin/env python
# coding=utf-8
import sys
import gzip
from Bio import SeqIO

Genome = gzip.open(sys.argv[1], 'rt')
minLength = int(sys.argv[2])
Output = gzip.open("sup2kb" + sys.argv[1], 'wt')

N = 0
big_scaffolds = []

print("Reading genome...")
for scaffold in SeqIO.parse(Genome, "fasta"):
    N += 1
    if len(scaffold.seq) >= minLength:
        big_scaffolds.append(scaffold)

print("Found", len(big_scaffolds), "scaffolds larger than", minLength, "bp within a total of", N)
print("Writing output...")
SeqIO.write(big_scaffolds, Output, "fasta")

Genome.close()
Output.close()
