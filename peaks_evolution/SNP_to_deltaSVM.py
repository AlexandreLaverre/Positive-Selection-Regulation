#!/usr/bin/env python
# coding=utf-8

import sys
from Bio import SeqIO
import pandas as pd
import gzip

#path = '/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/'
#vcf_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/VCF/filtered_chr21.vcf.gz'
#deltaSVM_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/deltas/ancestral_all_possible_deltaSVM.txt'
#ancestral_seq_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/sequences/filtered_ancestral_sequences.fa'
#maxLL_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/MLE_summary_50bins.csv'
#output_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/SNP_to_deltaSVM.txt'

VCF = pd.read_csv(sys.argv[1], sep='\t', header=None, compression='gzip', skiprows=20)
DeltaSVM = pd.read_csv(sys.argv[2], sep='\t', header=0)
AncestralSeq = SeqIO.to_dict(SeqIO.parse(open(sys.argv[3]), "fasta"))
MaxLL = pd.read_csv(sys.argv[4], header=0)
output = open(sys.argv[5], 'w')

# Correctly retrieve the header
with gzip.open(sys.argv[1], 'rt') as file:
    for line in file:
        if line.startswith('#CHROM'):
            header = line.strip('\t').split('\t')
            break
VCF.columns = header + ["chr", "start", "end", "PeakID"]


# Write header
tab = '\t'
pop = ['EAS_AF', 'EUR_AF', 'AFR_AF', 'AMR_AF']
output.write("ID\tPos\tRef\tAlt\t"+'\t'.join(pop)+"\tDeltaSVM\tStabParam\tAlphaPos\tBetaPos\n")

for SNP in VCF.iterrows():
    ID = SNP[1]['PeakID']
    ref = SNP[1]['REF']
    alt = SNP[1]['ALT']
    freq = [pop.split('=')[1] for pop in SNP[1]['INFO'].split(';')[4:8]]

    # Relative position in sequence
    start, end = int(ID.split('_')[0].split(':')[1]), int(ID.split('_')[0].split(':')[2])
    pos = int(SNP[1]['POS']) - start

    # Check that an alignment exists
    if ID not in AncestralSeq.keys():
        continue

    # Check that sequence do not contain gaps
    original_length = end-start
    aligned_length = len(AncestralSeq[ID].seq)

    if aligned_length != original_length:
        continue

    # Check that MaxLL estimations exist for this sequence
    if ID not in MaxLL['ID'].values:
        continue
    else:
        stab_param = MaxLL.loc[MaxLL['ID'] == ID, "AlphaPurif"].iloc[0]
        pos_params = MaxLL.loc[MaxLL['ID'] == ID, ["AlphaPos", "BetaPos"]].iloc[0]

    # Check that ref or alt allele is in ancestral sequence and get deltaSVM
    allSVM = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0]

    # Ref = Ancestral
    if pd.isna(allSVM[f"pos{pos}:{ref}"]):
        SNP_deltaSVM = allSVM[f"pos{pos}:{alt}"]
    # Alt = Ancestral
    elif pd.isna(allSVM[f"pos{pos}:{alt}"]):
        SNP_deltaSVM = -allSVM[f"pos{pos}:{ref}"]
        freq = [str(1-float(pop)) for pop in freq]
    else:
        continue

    # Write output
    output.write(f"{ID}\t{pos}\t{ref}\t{alt}\t{tab.join(freq)}\t{SNP_deltaSVM}\t{stab_param}"
                 f"\t{pos_params['AlphaPos']}\t{pos_params['BetaPos']}\n")
