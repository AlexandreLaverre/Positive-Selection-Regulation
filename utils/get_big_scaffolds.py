#!/usr/bin/env python
# coding=utf-8
import sys
import gzip
from Bio import SeqIO

Genome = gzip.open(sys.argv[1], 'rt')
minLength = int(sys.argv[2])
Output = gzip.open("sup2kb_" + sys.argv[1], 'wt')
names = open("sup2kb_" + sys.argv[1] + "_scaffold_names.txt", 'wt')
N = 0
GenomeSize = 0
big_scaffolds = []

print("Reading genome...")
for scaffold in SeqIO.parse(Genome, "fasta"):
    N += 1
    GenomeSize += len(scaffold)
    if len(scaffold.seq) >= minLength:
        big_scaffolds.append(scaffold)
        names.write(scaffold.id + '\n')

Genome.close()
names.close()

GenomeSize = round(GenomeSize/1000000000, 2)
FilteredGenomeSize = round(sum([len(scaffold) for scaffold in big_scaffolds])/1000000000, 2)
PropFiltered = round((FilteredGenomeSize/GenomeSize)*100, 2)

print("Found", len(big_scaffolds), "scaffolds larger than", minLength, "bp within a total of", N)
print("Filtered genome size = ", FilteredGenomeSize, "Gb over", GenomeSize, "Gb (", PropFiltered, "%)")
print("Writing and zipping output...")
SeqIO.write(big_scaffolds, Output, "fasta")

Output.close()
