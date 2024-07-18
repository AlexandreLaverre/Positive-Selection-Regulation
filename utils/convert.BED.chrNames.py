#!/usr/bin/env python
# coding=utf-8
import sys
import os

sp = sys.argv[1]
sample = sys.argv[2]
TF = sys.argv[3]

path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel" if sys.argv[4] == "cluster" \
    else "/Users/alaverre/Documents/Detecting_positive_selection/"

Correspondence = f"{path}/data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
pathPeaks = f"{path}/results/peaks_calling/NarrowPeaks/{sp}/{sample}/"
pathMatrix = f"{path}/results/substitution_matrix/{sp}/"
chromosomes = [os.path.splitext(file)[0] for file in os.listdir(pathMatrix)]

####################################################################################################
Correspondence_dict = {}
with open(Correspondence, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        Correspondence_dict[i[0]] = i[1]

for data in ["peaks", "consensus_summits"]:
    BED = f"{pathPeaks}/{TF}.{data}.bed"
    if not os.path.exists(BED):
        continue
    outfile = open(f"{pathPeaks}/{TF}.{data}_UCSC_names.bed", 'w')
    with open(BED, 'r') as f2:
        for i in f2.readlines():
            i = i.strip("\n")
            i = i.split("\t")

            old = str(i[0])
            peak_ID = str(i[3]) if data == "peaks" else str(i[3].split(":")[3])
            if old in Correspondence_dict.keys():
                new_chr = str(Correspondence_dict[old])

                # Filter chromosomes based on substitution matrix.
                if old in chromosomes or new_chr in chromosomes:
                    new_ID = f"{new_chr}:{str(i[1])}:{str(i[2])}_{peak_ID}"
                    outfile.write("\t".join([new_chr, str(i[1]), str(i[2]), new_ID]) + '\n')

    outfile.close()
####################################################################################################
