#!/usr/bin/env python
# coding=utf-8

import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
import gzip
import Positive_Selection_Tests.functions.MLEvol as ML

path = '/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/'

vcf_file = sys.argv[1] if len(sys.argv) > 1 else path + 'polymorphism_analyses/NarrowPeaks/drosophila/modERN/CG13204_CG13204-GFP_pupa_1/filtered_chrX.vcf.gz'
deltaSVM_file = sys.argv[2] if len(sys.argv) > 2 else path + 'positive_selection/NarrowPeaks/drosophila/modERN/CG13204_CG13204-GFP_pupa_1/deltas/focal_ancestral_all_possible_deltaSVM.txt'
focal_seq_file = sys.argv[3] if len(sys.argv) > 3 else path + 'positive_selection/NarrowPeaks/drosophila/modERN/CG13204_CG13204-GFP_pupa_1/sequences/filtered_focal_ancestral_sequences.fa'
genome_file = sys.argv[4] if len(sys.argv) > 4 else path + "../data/genome_sequences/drosophila/dm6.fa.gz"
maxLL_file = sys.argv[5] if len(sys.argv) > 5 else path + 'positive_selection/NarrowPeaks/drosophila/modERN/CG13204_CG13204-GFP_pupa_1/Tests/MLE_summary_exact_ranked_ancestral.csv'
output_file = sys.argv[6] if len(sys.argv) > 6 else path + 'polymorphism_analyses/NarrowPeaks/drosophila/modERN/CG13204_CG13204-GFP_pupa_1/SNP_to_deltaSVM_chrX.txt'

print("Reading input files...")
DeltaSVM = pd.read_csv(deltaSVM_file, sep='\t', header=0)
#FocalSeq = SeqIO.to_dict(SeqIO.parse(open(focal_seq_file), "fasta"))
genome = SeqIO.to_dict(SeqIO.parse(gzip.open(genome_file, "rt"), "fasta"))
MaxLL = pd.read_csv(maxLL_file, header=0)
output = open(output_file, 'w')
sp = genome_file.split("/")[-2]

# Correctly retrieve VCF header
with gzip.open(vcf_file, 'rt') as file:
    for i, line in enumerate(file):
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            row_to_skip = i
            break

header.extend(["chr", "start", "end", "PeakID"])
VCF = pd.read_csv(vcf_file, sep='\t', header=None, names=header, compression='gzip', skiprows=row_to_skip+1)

# Write header
tab = '\t'
pop = ['EAS_AF', 'EUR_AF', 'AFR_AF', 'AMR_AF', 'SAS_AF'] if sp == "human" else ["DGRP2"]
output.write("ID\tPos\tRef\tAlt\tNbAlt\tNbTot\t" + '\t'.join(pop) + "\tLength\tFlag\t" +
             "DeltaSVM\tSelCoefStab\tProbabStab\tSelCoefPos\tProbabPos\n")


print("Analysing VCF file...")
tot, valid, noMax, noDelta, indel = 0, 0, 0, 0, 0

for idx, SNP in VCF.iterrows():
    tot += 1
    print(sp, SNP['ID'])
    # Skip indels
    if sp != "human" and "SNP" not in SNP['ID']:
        indel += 1
        continue

    ID = SNP['PeakID']
    print(tot, ID)
    ID_delta = ID if sp == "drosophila" else f"{ID.split('_')[0]}:{vcf_file.split('/')[-2]}"

    # Check that focal sequence exists for this ID
    #if ID not in FocalSeq.keys():
    #   continue

    chr = SNP['#CHROM']
    SNP_pos = int(SNP['POS'] - 1)
    ref, alt = SNP['REF'], SNP['ALT']
    assert genome[chr].seq[SNP_pos].upper() == ref, "Reference does not correspond to the genome sequence."

    # Get allele frequencies
    nb_alt, nb_tot = SNP['INFO'].split(';')[0].split('=')[1], SNP['INFO'].split(';')[1].split('=')[1]
    if sp == "human":
        freq = [pop.split('=')[1] for pop in SNP['INFO'].split(';')[4:9]]  # frequencies for each population
    else:
        nb_tot = str(int(nb_tot)+int(nb_alt))  # in drosophila nb_tot correspond to nb of reference alleles
        freq = [str(int(nb_alt) / int(nb_tot))]

    # Filter out sequences that are too short or too long
    start, end = int(ID.split('_')[0].split(':')[1]), int(ID.split('_')[0].split(':')[2])
    length = end - start
    if length < 20 or length > 1000:
        continue

    # Check that MaxLL estimations exist for this sequence
    if ID not in MaxLL['ID'].values:
        noMax += 1
        continue
    else:
        Stab_param = MaxLL.loc[MaxLL['ID'] == ID, "AlphaPurif"].iloc[0]
        Pos_params = MaxLL.loc[MaxLL['ID'] == ID, ["AlphaPos", "BetaPos"]].iloc[0]

    # Get deltaSVM calculated from Ancestral sequence
    if ID_delta not in DeltaSVM['ID'].values:
        noDelta += 1
        continue

    allSVM = DeltaSVM.loc[DeltaSVM['ID'] == ID_delta, "pos0:A":].iloc[0]
    allSVM_noNA = allSVM.dropna().values.tolist()

    seq_ids = allSVM[allSVM.isna()].index.tolist()[0:length]
    allSVM_ID = ML.get_svm_ids(seq_ids)

    sorted_SVM, x = zip(*sorted(zip(allSVM_noNA, allSVM_ID)))

    # Get scaled bins reusing MaxLL function
    y, z, scaled_bins = ML.get_svm_exact(allSVM_noNA, allSVM_noNA[1], allSVM_ID, "NA", norm="ranked", get_mut_rate=False)

    # Relative position in sequence
    pos = SNP_pos - start
    # Get ID from posSet instead of focal_sequences !!
    #assert genome[chr].seq[SNP_pos].upper() == FocalSeq[ID].seq[pos], "SNP position does not correspond to the genome sequence."

    # Ref = Ancestral
    if pd.isna(allSVM[f"pos{pos}:{ref}"]):
        SNP_deltaSVM = allSVM[f"pos{pos}:{alt}"]
        flag = "ref_ancestral"
    # Alt = Ancestral
    elif pd.isna(allSVM[f"pos{pos}:{alt}"]):
        SNP_deltaSVM = -allSVM[f"pos{pos}:{ref}"]
        flag = "alt_ancestral"
    else:
        #print(f"Weird! SNP not found in ALlSVM. Ref:{ref}; Alt: {alt}; sequence: {FocalSeq[ID].seq[pos]}")
        continue

    # Get SNP index
    SNP_index = ML.get_obs_index(SNP_deltaSVM, sorted_SVM, flag)[0]

    # Compute selection coefficient
    Coef_Stab = ML.coeff_selection(scaled_bins[SNP_index], np.array([Stab_param]))
    Coef_Pos = ML.coeff_selection(scaled_bins[SNP_index], np.array([Pos_params['AlphaPos'], Pos_params['BetaPos']]))

    # Compute probabilities of fixation
    Probab_Stab = ML.proba_fixation(Coef_Stab)
    Probab_Pos = ML.proba_fixation(Coef_Pos)

    # Write output
    valid += 1
    output.write(f"{ID}\t{pos}\t{ref}\t{alt}\t{nb_alt}\t{nb_tot}\t{tab.join(freq)}\t{length}\t{flag}\t"
                 f"{SNP_deltaSVM}\t{Coef_Stab}\t{Probab_Stab}\t{Coef_Pos}\t{Probab_Pos}\n")


output.close()

print(f"Total SNPs: {tot}; Valid SNPs: {valid}; Indel: {indel}; No Max estimation: {noMax}; noDelta: {noDelta}.\n")