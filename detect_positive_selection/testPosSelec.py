#!/usr/bin/env python
# coding=utf-8

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas
from alive_progress import alive_bar

species = sys.argv[1]       # Species name: human dog
NbRand = int(sys.argv[2])        # Number of random substitutions permutations per sequence
Evol = sys.argv[3]          # Substitution model: uniform or proba

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
pathSelection = path + "positive_selection/" + species + "/CEBPA/"
Ancestral_fasta = pathSelection + "filtered_ancestral_sequences.fa_subsample"
Focal_fasta = pathSelection + "filtered_focal_sequences.fa_subsample"
ModelEstimation = pathSelection + "Model/kmer_predicted_weight.txt"
pathSubMat = path + "/substitution_matrix/dog/"

Output = open(pathSelection + "PosSelTest_deltaSVM_" + str(NbRand) + "permutations.txt_subsample_new", "w")


####################################################################################################
# Get number of substitutions per sequence
def get_sub_number(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq1, seq2))


# Calculate delta SVM from sliding windows
def calculate_delta_svm(seq_ref, seq_alt, kmer_len, svm_scores):
    delta_svm = 0
    for pos in range(len(seq_ref) - kmer_len + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos+kmer_len]
        kmer_alt = seq_alt[pos:pos+kmer_len]
        delta_svm += svm_scores[kmer_alt] - svm_scores[kmer_ref]  # sum of delta between sequences for each kmer
    return delta_svm


# Get random sequences according to substitution matrix
def get_random_seqs(seq, sub_prob, sub_prob_norm, nb_sub, nb_perm):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nucl = seq[pos]                                # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nucl].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of all positions to sum at 1

    random_seqs = []
    for perm in range(nb_perm):
        rand_seq = seq
        # draw positions
        rand_pos = np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, replace=False, size=nb_sub)

        # draw directions
        for pos in rand_pos:
            old_nuc = seq[pos]
            directions = list(sub_prob_norm[old_nuc].keys())
            proba = list(sub_prob_norm[old_nuc].values())
            new_nuc = np.random.choice(directions, p=proba)
            rand_seq = rand_seq[:pos] + new_nuc + rand_seq[pos+1:]  # change value in pos

        random_seqs.append(rand_seq)

    return random_seqs


####################################################################################################
# Binding affinity values per kmer
print("Get SVM score for each kmer from:", ModelEstimation)
SVM_dict = {}
with open(ModelEstimation, 'r') as model:
    for i in model.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        kmer = i[0]
        KmerLen = len(kmer)
        svm_score = float(i[1])
        rev_kmer = str(Seq(kmer).reverse_complement())  # add the reverse complement kmer

        SVM_dict[kmer] = svm_score
        SVM_dict[rev_kmer] = svm_score


# Substitution matrix
if Evol == 'uniform':
    SubMat_uniform = {'A': {'A': 0, 'C': 0.333, 'G': 0.333, 'T': 0.334},
                      'C': {'A': 0.333, 'C': 0, 'G': 0.333, 'T': 0.334},
                      'G': {'A': 0.333, 'C': 0.333, 'G': 0, 'T': 0.334},
                      'T': {'A': 0.333, 'C': 0.333, 'G': 0.334, 'T': 0}}
else:
    print("Get substitution matrix for each chromosome...")
    SubMats = {}
    SubMats_norm = {}
    for file in os.listdir(pathSubMat):
        if file.endswith('.txt'):
            chrom = 'chr' + file.strip('.txt')

            chrom_Table = pandas.read_table(pathSubMat + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']     # change row values
            np.fill_diagonal(chrom_Table.values, 0)      # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            SubMats[chrom] = chrom_SubMat

            # Get substitution probabilities of each nucleotide to sum at 1 for drawing directions
            chrom_Table_norm = chrom_Table.div(chrom_Table.sum(axis=1), axis=0)
            chrom_SubMat_norm = chrom_Table_norm.to_dict('index')
            SubMats_norm[chrom] = chrom_SubMat_norm


# Sequences
print("Get Ancestral sequences...")
AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))

print("Get Focal sequences...")
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))

####################################################################################################
Output.write("ID\tdeltaSVM\tNbSub\tpval.high\n")  # header

print("Running Positive Selection test for each sequence with", NbRand, "random permutations...")
with alive_bar(len(FocalSeqs.keys())) as bar:
    for ID in FocalSeqs.keys():
        bar()   # Progress bar
        FocalSeq = str(FocalSeqs[ID].seq)
        AncestralSeq = str(AncestralSeqs[ID].seq)

        # Get corresponding substitution matrix
        chrom = ID.split(':')[0]
        SubMat_proba = SubMats[chrom] if Evol != 'uniform' else SubMat_uniform
        SubMat_proba_normed = SubMats_norm[chrom] if Evol != 'uniform' else SubMat_uniform

        # Number of substitutions between Ancestral and Focal sequences
        NbSub = get_sub_number(FocalSeq, AncestralSeq)
        if NbSub > 1:

            DeltaObs = calculate_delta_svm(AncestralSeq, FocalSeq, KmerLen, SVM_dict)
            RandomSeqs = get_random_seqs(AncestralSeq, SubMat_proba, SubMat_proba_normed, NbSub, NbRand)
            DeltaRand = [calculate_delta_svm(AncestralSeq, RandSeq, KmerLen, SVM_dict) for RandSeq in RandomSeqs]

            # Calculate p-value
            NbHigherRand = sum(rand > DeltaObs for rand in DeltaRand)
            pvalHigh = NbHigherRand / len(DeltaRand)

            Output.write(ID + "\t" + str(DeltaObs) + "\t" + str(NbSub) + "\t" + str(pvalHigh) + "\n")

Output.close()

print("Done!")

####################################################################################################
