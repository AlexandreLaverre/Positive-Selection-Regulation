#!/usr/bin/env python
# coding=utf-8
import sys
import gzip

species = sys.argv[1]
GFF = sys.argv[2]
suffix = sys.argv[3]
cluster = sys.argv[4]
reverse = sys.argv[5]

if cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/data/genome_sequences/" + species
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/data/genome_sequences/" + species

Correspondence = path + "/chromosome_correspondence_" + suffix + ".txt"
output = gzip.open(f"{path}/{suffix}_{GFF}", 'wt')

####################################################################################################
Correspondence_dict = {}
with open(Correspondence, 'r') as f1:
    for i in f1.readlines():
        i = i.split("\t")
        old = i[0] if reverse == "F" else i[1]
        new = i[1] if reverse == "F" else i[0]
        Correspondence_dict[old] = new

missing_names = []
with gzip.open(f"{path}/{GFF}", 'rt') as f2:
    for line in f2:
        if line.startswith('#'):
            output.write(line)
        else:
            line = line.split("\t")
            old_name = str(line[0])
            try:
                new_name = str(Correspondence_dict[old_name])
                output.write(new_name + '\t' + '\t'.join(line[1:]))
            except:
                missing_names.append(old_name)
                continue

print("Found", len(set(missing_names)), "missing scaffolds.")
print("First 10 missing:"', '.join(list(set(missing_names))[0:10]))
print("Done!")

output.close()
####################################################################################################
