#!/usr/bin/env python
# coding=utf-8

from Bio import SeqIO
import numpy as np
from sklearn.preprocessing import normalize

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/human/CEBPA/"
Ancestral_fasta = path + "filtered_ancestral_sequences.fa"
Focal_fasta = path + "filtered_focal_sequences.fa"
ModelEstimation = path + "Model/kmer_predicted_weight.txt"
Substitution_file = path
NbRand = 10
KmerLen = 11
Output = open(path + "PosSelTest_deltaSVM_" + str(NbRand) + "permutations.txt")


####################################################################################################
# Get number of substitutions per sequence
def get_sub_number(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq1, seq2))


def calculate_delta_svm(seq1, seq2, kmer_len, svm_scores):
    delta_svm = pos = 0
    while pos < len(seq1) - kmer_len:   # sliding window of kmer length
        kmer_ref = seq1[pos:kmer_len]
        kmer_alt = seq2[pos:kmer_len]
        delta_svm += svm_scores[kmer_ref] - svm_scores[kmer_alt]  # sum of delta between sequences for each kmer
        pos += 1
    return delta_svm


# Get random sequences from substitution matrix
def get_random_seqs(seq, mat_prob, mat_prob_norm, mat_dir, nb_sub, nb_perm):
    # Get substitution probability for each position:
    # the probability of a substitution at a given position is the sum
    # of the probabilities to get a substitution from this nucleotide
    pos_proba = np.empty(shape=(len(seq), 1), dtype=float)
    for pos in range(len(seq)):
        nuc = seq[pos]
        pos_proba[pos] = sum(mat_prob[nuc])
    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of positions to sum at 1

    random_seqs = []
    for perm in range(1, nb_perm):
        rand_seq = seq
        rand_pos = np.random.choice(normed_pos_proba.size, p=normed_pos_proba, replace=False, size=nb_sub)  # draw position
        for pos in rand_pos:
            old_nuc = seq[pos]
            new_nuc = np.random.choice(mat_dir[old_nuc], p=mat_prob_norm[old_nuc], size=1)  # draw direction
            rand_seq = rand_seq[:rand_pos] + new_nuc + rand_seq[rand_pos+1:]  # change value in pos

        random_seqs.append(rand_seq)

    return random_seqs


# Get p-values from distribution of random delta SVM
def test_delta_svm(delta_obs, delta_rand):
    nb_higher = nb_lower = 0
    for rand in delta_rand:
        if delta_obs >= rand: nb_higher += 1
        if delta_obs <= rand: nb_lower += 1

    pval_higher = nb_higher/len(delta_rand)
    pval_lower = nb_lower/len(delta_rand)

    return pval_higher, pval_lower


####################################################################################################
# Get SVM score for each kmer
SVM_dict = {}
with open(ModelEstimation, 'r') as model:
    for i in model.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        kmer = i[0]
        svm_score = i[1]
        rev_kmer = kmer.reverse_complement()

        SVM_dict[kmer] = svm_score
        SVM_dict[rev_kmer] = svm_score

# Get substitutions matrix
#SubMat = open(Substitution_file)
SubMat_direction = {'A': ['T', 'C', 'G'], 'T': ['A', 'C', 'G'], 'C': ['A', 'T', 'G'], 'G': ['A', 'T', 'C']}
SubMat_proba = {'A': [0.7, 0.15, 0.15], 'T': [0.7, 0.15, 0.15], 'C': [0.7, 0.15, 0.15], 'G': [0.7, 0.15, 0.15]}
SubMat_proba_normed = {}
for nuc in SubMat_proba.keys():
    SubMat_proba_normed[nuc] = SubMat_proba[nuc] / sum(SubMat_proba[nuc])  # normalisation of direction to sum at 1

# Get sequences for each ID
AncestralSeq = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))
FocalSeq = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))

####################################################################################################
# Running test for each sequence
Output.write("ID\tdeltaSVM\tNbSub\tpval.low\tpval.high")  # header

for ID in FocalSeq.keys():
    NbSub = get_sub_number(FocalSeq[ID], AncestralSeq[ID])

    if NbSub > 1:
        DeltaObs = calculate_delta_svm(AncestralSeq[ID], FocalSeq[ID], KmerLen, SVM_dict)

        RandomSeqs = get_random_seqs(FocalSeq[ID], SubMat_proba, SubMat_proba_normed, SubMat_direction, NbSub, NbRand)
        DeltaRand = [calculate_delta_svm(AncestralSeq[ID], Seq, KmerLen, SVM_dict) for Seq in RandomSeqs]

        pvalHigh, pvalLow = test_delta_svm(DeltaObs, DeltaRand)

        Output.write(ID + str(DeltaObs) + str(NbSub) + str(pvalHigh) + str(pvalLow) + "\n")

####################################################################################################


# Problem: possible to draw the same position (with different direction) several times
# define proba of position first = sum(proba direction), then draw the direction
"""matrix = np.empty(shape=(len(seq), 3), dtype=float)
for pos in range(len(seq)):
    nuc = seq[pos] 
    matrix[pos] = mat_prob[nuc]  
normed_matrix = matrix / sum(sum(matrix))  # normalisation to sum at 1

random_seqs = []
# Draw substitutions according to probabilities
for perm in range(1, nb_perm):
    rand_seq = seq
    rand_sub = np.random.choice(np.arange(normed_matrix.size), p=normed_matrix.ravel(), replace=False, size=nb_sub)
"""

