#!/usr/bin/env python3
# coding=utf-8
import sys
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

###################################################################################################
path = "/Users/alaverre/Documents/Detecting_positive_selection/data/genome_alignments/Canidae/per_chrom/"
InputAlignFile = path + "FASTAs/" + sys.argv[1]
OutputAlignFile = path + "new_PHYLIPs/" + sys.argv[2]

###################################################################################################

align = AlignIO.read(InputAlignFile, "fasta")

ids = []
seqs = []
for seq in align:
    ids.append(seq.id)
    seqs.append(seq.seq)

records = (SeqRecord(s, id=i) for (i, s) in zip(ids, seqs))

with open(OutputAlignFile, "w") as OutFile:
    AlignIO.write(MultipleSeqAlignment(records), OutFile, "phylip-relaxed")


###################################################################################################
