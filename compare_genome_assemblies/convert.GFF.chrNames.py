#!/usr/bin/env python
# coding=utf-8
import sys
import gzip

species = sys.argv[1]
GFF = sys.argv[2]
suffix = sys.argv[3]
cluster = sys.argv[4]

if cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/"

Correspondence = path + "data/genome_sequences/" + species + "/chromosome_correspondence_" + suffix + ".txt"

output = open(GFF + "_UCSC_names", 'w')

####################################################################################################
Correspondence_dict = {}
with open(Correspondence, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        Correspondence_dict[i[0]] = i[1]

with gzip.open(GFF, 'r') as f2:
    for line in f2:
        if line.startswith('#'):
            output.write(line)
        else:
            line = line.split("\t")

            old_name = str(line[0])
            new_name = str(Correspondence_dict[old_name])
            output.write(new_name + '\t' + '\t'.join(i[1:]))

output.close()
####################################################################################################
