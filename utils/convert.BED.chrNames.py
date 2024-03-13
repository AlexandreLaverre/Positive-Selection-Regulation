#!/usr/bin/env python
# coding=utf-8
import sys
import os

path = os.getcwd() + "/"

BED = path + sys.argv[1]
Correspondence = path + sys.argv[2]
output = path + sys.argv[3]

outfile = open(output, 'w')
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

        old = str(i[0])
        if old in Correspondence_dict.keys():
            new_chr = str(Correspondence_dict[old])
            new_ID = new_chr + ':' + str(i[1]) + ':' + str(i[2])

            outfile.write(new_chr + '\t' + str(i[1]) + '\t' + str(i[2]) + '\t' + new_ID + '\n')

outfile.close()
####################################################################################################
