#!/usr/bin/env python
# coding=utf-8

import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
import gzip
import Positive_Selection_Tests.functions.MLEvol as ML

# path = '/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/'
# vcf_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/VCF/filtered_chr21.vcf.gz'
# deltaSVM_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/deltas/ancestral_all_possible_deltaSVM.txt'
# focal_seq_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/sequences/filtered_focal_sequences.fa'
# genome_file = "/Users/alaverre/Documents/Detecting_positive_selection/data/genome_sequences/human/hg38.fa.gz"
# maxLL_file = path + 'positive_selection/NarrowPeaks/human/Wilson/CEBPA/MLE_summary_50bins.csv'
# output_file = path + 'polymorphism_analyses/NarrowPeaks/human/Wilson/CEBPA/SNP_to_deltaSVM_chr21.txt'

# Read input files
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
output.write("ID\tPos\tRef\tAlt\tNbAlt\tNbTot\t" + '\t'.join(pop) + "\tLength\tFlag\t" +
             "DeltaSVM\tConclusion\tStabParam\tSelCoefStab\tProbabStab\t" +
             "AlphaPos\tBetaPos\tSelCoefPos\tProbabPos\n")

for SNP in VCF.iterrows():
    ID = SNP[1]['PeakID']
    chr = SNP[1]['#CHROM']
    SNP_pos = int(SNP[1]['POS'] - 1)
    ref, alt = SNP[1]['REF'], SNP[1]['ALT']
    freq = [pop.split('=')[1] for pop in SNP[1]['INFO'].split(';')[4:8]]
    nb_alt, nb_tot = SNP[1]['INFO'].split(';')[0].split('=')[1], SNP[1]['INFO'].split(';')[1].split('=')[1]

    assert genome[chr].seq[SNP_pos].upper() == ref, "Reference does not correspond to the genome sequence."

    # Check that focal sequence exists for this ID
    if ID not in FocalSeq.keys():
        continue

    # Check that sequence do not contain gaps
    start, end = int(ID.split('_')[0].split(':')[1]), int(ID.split('_')[0].split(':')[2])
    original_length = end - start
    aligned_length = len(FocalSeq[ID].seq)
    if aligned_length != original_length:
        continue

    assert genome[chr].seq[start].upper() == FocalSeq[ID].seq[0], "Start does not correspond to the genome sequence."

    # Check that MaxLL estimations exist for this sequence
    if ID not in MaxLL['ID'].values:
        continue
    else:
        Stab_param = MaxLL.loc[MaxLL['ID'] == ID, "AlphaPurif"].iloc[0]
        Pos_params = MaxLL.loc[MaxLL['ID'] == ID, ["AlphaPos", "BetaPos"]].iloc[0]
        Conclusion = MaxLL.loc[MaxLL['ID'] == ID, "Conclusion"].iloc[0]

    # Relative position in sequence
    pos = SNP_pos - start
    assert genome[chr].seq[SNP_pos].upper() == FocalSeq[ID].seq[
        pos], "Position does not correspond to the genome sequence."

    # Get deltaSVM calculated from Ancestral sequence
    allSVM = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0]

    # Ref = Ancestral
    if pd.isna(allSVM[f"pos{pos}:{ref}"]):
        SNP_deltaSVM = allSVM[f"pos{pos}:{alt}"]
        flag = "ref_ancestral"
    # Alt = Ancestral
    elif pd.isna(allSVM[f"pos{pos}:{alt}"]):
        SNP_deltaSVM = -allSVM[f"pos{pos}:{ref}"]
        flag = "alt_ancestral"
    else:
        continue

    # Compute selection coefficient
    SVM_bounds = [np.nanmin(allSVM), np.nanmax(allSVM)]
    Coef_Stab = ML.coeff_selection_s(SNP_deltaSVM, [Stab_param], SVM_bounds)
    Coef_Pos = ML.coeff_selection_s(SNP_deltaSVM, [Pos_params['AlphaPos'], Pos_params['BetaPos']], SVM_bounds)
    Probab_Stab = ML.proba_fixation_s(Coef_Stab)
    Probab_Pos = ML.proba_fixation_s(Coef_Pos)

    # Write output
    output.write(f"{ID}\t{pos}\t{ref}\t{alt}\t{nb_alt}\t{nb_tot}\t{tab.join(freq)}\t{original_length}\t{flag}\t"
                 f"{SNP_deltaSVM}\t{Conclusion}\t{Stab_param}\t{Coef_Stab}\t{Probab_Stab}\t"
                 f"{Pos_params['AlphaPos']}\t{Pos_params['BetaPos']}\t{Coef_Pos}\t{Probab_Pos}\n")

output.close()
