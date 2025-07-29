#!/usr/bin/env python
# coding=utf-8
import sys
import os

cwd = os.getcwd()

pathCorrespondence = f"{cwd}/{sys.argv[1]}"
pathBED = f"{cwd}/{sys.argv[2]}"
output = open(f"{cwd}/{sys.argv[3]}", 'w')

Correspondence_dict = {}
with open(pathCorrespondence, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        Correspondence_dict[i[0]] = i[1]

with open(pathBED, 'r') as f1:
    for i in f1.readlines():
        i = i.strip("\n").split("\t")
        old = str(i[0])

        # skip header
        if old == "chrom":
            continue

        if old in Correspondence_dict.keys():
            new_chr = str(Correspondence_dict[old])

        output.write("\t".join([new_chr, str(i[1]), str(i[2]), str(i[3])]) + '\n')

output.close()
####################################################################################################
