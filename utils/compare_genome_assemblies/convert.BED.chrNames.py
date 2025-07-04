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
chromosomes_list = [chrom.strip('.txt') for chrom in os.listdir(MatrixPath)]

if not os.path.isfile(Correspondence):
    sys.exit("Chromosome correspondence file doesn't exist!")

output = open(BED + "_UCSC_names", 'w')
output_interval = open(BED + "_UCSC_names_interval", 'w')

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
        interval = i[3].split(":")[3]

        # Remove ID in scaffolds
        old = str(i[0])
        old_chr = 'chr' + str(i[0])
        if old in Correspondence_dict.keys():
            new_chr = str(Correspondence_dict[old])

            # get only the main chromosome and not all the scaffolds
            if any(chr in chromosomes_list for chr in [old, old_chr, new_chr]):
                new_ID = new_chr + ':' + str(i[1]) + ':' + str(i[2])

                output.write(new_chr + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + new_ID + '\n')
                output_interval.write(new_chr + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + new_ID + '\t' + interval + '\n')

output.close()
####################################################################################################
