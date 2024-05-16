#!/usr/bin/env python
# coding=utf-8

import sys
from Bio import SeqIO
import pandas as pd
import gzip

#path = '/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/'
#vcf_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/VCF/filtered_chr21.vcf.gz'
#deltaSVM_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/deltas/ancestral_all_possible_deltaSVM.txt'
#seq_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/sequences/filtered_focal_sequences_upper.fa'
#maxLL_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/MLE_summary_50bins.csv'
#output_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/SNP_to_deltaSVM.txt'
#genome_file = "/Users/alaverre/Documents/Detecting_positive_selection/data/genome_sequences/human/hg38.fa.gz"

VCF = pd.read_csv(sys.argv[1], sep='\t', header=None, compression='gzip', skiprows=20)
DeltaSVM = pd.read_csv(sys.argv[2], sep='\t', header=0)
FocalSeq = SeqIO.to_dict(SeqIO.parse(open(sys.argv[3]), "fasta"))
genome = SeqIO.to_dict(SeqIO.parse(gzip.open(sys.argv[4], "rt"), "fasta"))
MaxLL = pd.read_csv(sys.argv[5], header=0)
output = open(sys.argv[6], 'w')

# Correctly retrieve VCF header
with gzip.open(sys.argv[1], 'rt') as file:
    for line in file:
        if line.startswith('#CHROM'):
            header = line.strip('\t').split('\t')
            break
VCF.columns = header + ["chr", "start", "end", "PeakID"]

# Write header
tab = '\t'
pop = ['EAS_AF', 'EUR_AF', 'AFR_AF', 'AMR_AF']
output.write("ID\tPos\tRef\tAlt\t"+'\t'.join(pop)+"\tDeltaSVM\tStabParam\tAlphaPos\tBetaPos\tLength\n")

for SNP in VCF.iterrows():
    ID = SNP[1]['PeakID']
    chr = SNP[1]['#CHROM']
    SNP_pos = int(SNP[1]['POS']-1)
    ref, alt = SNP[1]['REF'], SNP[1]['ALT']

    assert genome[chr].seq[SNP_pos].upper() == ref, "Reference does not correspond to the genome sequence."

    # Check that focal sequence exists for this ID
    if ID not in FocalSeq.keys():
        continue

    # Check that sequence do not contain gaps
    start, end = int(ID.split('_')[0].split(':')[1]), int(ID.split('_')[0].split(':')[2])
    original_length = end-start
    aligned_length = len(FocalSeq[ID].seq)
    if aligned_length != original_length:
        continue

    assert genome[chr].seq[start].upper() == FocalSeq[ID].seq[0], "Start does not correspond to the genome sequence."

    # Check that MaxLL estimations exist for this sequence
    if ID not in MaxLL['ID'].values:
        continue
    else:
        stab_param = MaxLL.loc[MaxLL['ID'] == ID, "AlphaPurif"].iloc[0]
        pos_params = MaxLL.loc[MaxLL['ID'] == ID, ["AlphaPos", "BetaPos"]].iloc[0]

    # Relative position in sequence and frequency in populations
    pos = SNP_pos - start
    freq = [pop.split('=')[1] for pop in SNP[1]['INFO'].split(';')[4:8]]

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
                 f"\t{pos_params['AlphaPos']}\t{pos_params['BetaPos']}\t{original_length}\n")
