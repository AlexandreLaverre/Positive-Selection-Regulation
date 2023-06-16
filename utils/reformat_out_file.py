#!/usr/bin/env python

import sys
import gzip

Path_inFile = sys.argv[1]
Infile = gzip.open(Path_inFile, "rt") if Path_inFile.endswith(".gz") else open(Path_inFile, "r")

outFile = gzip.open(Path_inFile.strip(".gz") + "_2.gz", 'wt')

for line in Infile:
    clean_line = line.split()
    outFile.write("\t".join(clean_line) + '\n')

Infile.close()
outFile.close()

