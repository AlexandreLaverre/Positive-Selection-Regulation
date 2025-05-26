#!/usr/bin/env python
# coding=utf-8

import argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-sp", "--species", default="human", help="Species name: human or mouse")
parser.add_argument("-TF", "--sample", default="Wilson/HNF4A", help="Study name and Transcription Factor: Wilson/CEBPA Schmidt/CTCF ...")
parser.add_argument("-s", "--score", default="phastCons", help="phastCons or phyloP")
args = parser.parse_args()

path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/"
peaks = f"{path}/results/peaks_calling/NarrowPeaks/{args.species}/{args.sample}.peaks_UCSC_names.bed"
pathScore = f"{path}/results/{args.score}/NarrowPeaks/{args.species}/{args.sample}/"

if args.species == "mouse":
    suffix = f".{args.score}60way.glire.wigFix.gz.bed"
elif args.species == "drosophila":
    suffix = f"{args.score}27way.wigFix.bed.gz"
else:
    suffix = ".bed" if args.score == "phastCons" else ".phyloP17way.wigFix.gz.bed"

# Sort BED file
def sorted_dictionary(file):
    dic = defaultdict(list)
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip("\n").split("\t")
            chrom = line[0]
            if file == peaks:
                pos = (int(line[1]), int(line[2]), str(line[3]))
            else:
                pos = (int(line[1]), int(line[2]), str(line[4]))
            dic[chrom].append(pos)

    # Sorting coordinates for each chromosome
    if file == peaks:
        for k in dic.keys():
            dic[k] = list(set(tuple(x) for x in dic[k]))
            dic[k].sort(key=lambda x: x[0])

    return dic


peaks_dic = sorted_dictionary(peaks)

# Overlap interest to reference
dic_output = defaultdict(list)
count_ID_high = 0
count_score_high = 0
for chrom in peaks_dic.keys():
    score = f"{pathScore}/overlapping/{chrom}{suffix}"
    if not os.path.exists(score):
        # print(f"{chrom} doesn't exist in {args.score}!")
        continue
    else:
        print(chrom)

    score_dic = sorted_dictionary(score)
    first_i = 0
    for pos in peaks_dic[chrom]:
        start, end, ID = pos[0], pos[1], str(pos[2])
        # Initialization of first possible overlapping interest position
        i = first_i
        while i < len(score_dic[chrom]) and score_dic[chrom][i][1] < start:
            i += 1
        first_i = i

        # Adding all overlapping interest position to reference position
        while i < len(score_dic[chrom]) and score_dic[chrom][i][0] <= end:
            dic_output[ID].append(score_dic[chrom][i][2])
            i += 1

        # Check if lengths are equal
        ID_len = pos[1]-pos[0]
        score_len = len(dic_output[ID])
        if ID_len > score_len:
            count_ID_high += 1
        if ID_len < score_len:
            count_score_high += 1


print(f"Total peaks:{len(dic_output.keys())}; ID higher:{count_score_high}; Score higher:{count_ID_high}")

print("Writing output")
output = f"{path}/results/{args.score}/NarrowPeaks/{args.species}/{args.sample}/score_per_base.txt"
with open(output, 'w') as f:
    for ID in dic_output.keys():
        f.write(ID + "\t" + '\t'.join(dic_output[ID]) + "\n")


