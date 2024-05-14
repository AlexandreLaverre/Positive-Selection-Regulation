#!/usr/bin/env python
# coding=utf-8

import sys
from Bio import SeqIO
import pandas as pd

VCF = pd.read_csv(sys.argv[1], sep='\t', header=0)
DeltaSVM = pd.read_csv(sys.argv[2], sep='\t', header=0)
AncestralSeq = SeqIO.to_dict(SeqIO.parse(open(sys.argv[3]), "fasta"))
MaxLL = pd.read_csv(sys.argv[4], sep='\t', header=0)
output = open(sys.argv[5], 'w')

# Write header
output.write("ID\tPos\tRef\tAlt\tCount\tDeltaSVM\tStabParam\tAlphaPos\tBetaPos\n")

for SNP in VCF.iterrows():
    ID = SNP[1][-1]
    ref = SNP[1]['REF']
    alt = SNP[1]['ALT']
    count = SNP[1]['COUNT']

    # Relative position in sequence
    pos = SNP[1]['POS'] - ID.split(':')[1]

    # Check that sequence do not contain gaps
    original_length = ID.split(':')[2]-ID.split(':')[1]
    aligned_length = len(AncestralSeq[ID].seq)
    if aligned_length != original_length:
        continue

    # Check that MaxLL estimations exist for this sequence
    if ID not in MaxLL['ID'].values:
        continue
    else:
        stab_param = MaxLL.loc[MaxLL['ID'] == ID, "AlphaPurif"]
        pos_params = MaxLL.loc[MaxLL['ID'] == ID, ["AlphaPos", "BetaPos"]]

    # Check that ref or alt allele is in ancestral sequence and get deltaSVM
    allSVM = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0]

    # Ref = Ancestral
    if allSVM[f"pos{pos}:{ref}"] == "nan":
        SNP_deltaSVM = allSVM[f"pos{pos}:{alt}"]
    # Alt = Ancestral
    elif allSVM[f"pos{pos}:{alt}"] == "nan":
        SNP_deltaSVM = -allSVM[f"pos{pos}:{ref}"]
        count = 1-count
    else:
        continue

    # Write output
    output.write(f"{ID}\t{pos}\t{ref}\t{alt}\t{count}\t{SNP_deltaSVM}\t{stab_param}\t{pos_params[0]}\t{pos_params[1]}\n")
