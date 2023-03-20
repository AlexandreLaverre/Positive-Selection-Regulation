#!/usr/bin/env python
# coding=utf-8
import sys
import os

species = sys.argv[1]
BED = sys.argv[2]
cluster = sys.argv[3]

if cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/"

Correspondence = path + "data/genome_sequences/" + species + "/chromosome_correspondence.txt"
MatrixPath = path + "results/substitution_matrix/" + species
chromosomes_list = [chrom.split('.')[0] for chrom in os.listdir(MatrixPath)]

output = open(BED + "_UCSC_names", 'w')

####################################################################################################
Correspondence_dict = {}
with open(Correspondence, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        Correspondence_dict[i[0]] = i[1]

with open(BED, 'r') as f2:
    for i in f2.readlines():
        i = i.strip("\n")
        i = i.split("\t")

        # Remove ID in scaffolds
        old_chr = str(i[0])
        new_chr = str(Correspondence_dict[old_chr])
        if old_chr or new_chr in chromosomes_list:
            new_ID = new_chr + ':' + str(i[1]) + ':' + str(i[2])

            output.write(new_chr + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + new_ID + '\n')

output.close()
####################################################################################################
