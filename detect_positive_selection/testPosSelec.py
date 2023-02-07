#!/usr/bin/env python
# coding=utf-8

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from alive_progress import alive_bar

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/human/CEBPA/"
Ancestral_fasta = path + "filtered_ancestral_sequences.fa_subsample"
Focal_fasta = path + "filtered_focal_sequences.fa_subsample"
ModelEstimation = path + "Model/kmer_predicted_weight.txt"
Substitution_file = path  # unused for the moment
NbRand = 1000             # number of random substitutions permutations per sequence

Output = open(path + "PosSelTest_deltaSVM_" + str(NbRand) + "permutations.txt_subsample", "w")


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
def get_random_seqs(seq, sub_prob, sub_prob_norm, sub_dir, nb_sub, nb_perm):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nucl = seq[pos]                             # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nucl])        # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of all positions to sum at 1

    random_seqs = []
    for perm in range(nb_perm):
        rand_seq = seq

        # draw positions
        rand_pos = np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, replace=False, size=nb_sub)

        # draw directions
        for pos in rand_pos:
            old_nuc = seq[pos]
            new_nuc = np.random.choice(sub_dir[old_nuc], p=sub_prob_norm[old_nuc])
            rand_seq = rand_seq[:pos] + new_nuc + rand_seq[pos+1:]  # change value in pos

        random_seqs.append(rand_seq)

    return random_seqs


####################################################################################################
print("Get SVM score for each kmer from:", ModelEstimation)
SVM_dict = {}
with open(ModelEstimation, 'r') as model:
    for i in model.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        kmer = Seq(i[0])
        KmerLen = len(kmer)
        svm_score = float(i[1])
        rev_kmer = kmer.reverse_complement()  # add the reverse complement kmer

        SVM_dict[kmer] = svm_score
        SVM_dict[rev_kmer] = svm_score

# print("Get substitutions matrix from:", Substitution_file)
# SubMat = open(Substitution_file)
# Temporary uniform substitution matrix
SubMat_direction = {'A': ['T', 'C', 'G'], 'T': ['A', 'C', 'G'], 'C': ['A', 'T', 'G'], 'G': ['A', 'T', 'C']}
SubMat_proba = {'A': [0.33, 0.33, 0.33], 'T': [0.33, 0.33, 0.33], 'C': [0.33, 0.33, 0.33], 'G': [0.33, 0.33, 0.33]}

# Normalisation of probabilities to sum at 1 for each nucleotide --> useless if uniform
SubMat_proba_normed = {}
for nuc in SubMat_proba.keys():
    SubMat_proba_normed[nuc] = [sub_proba / sum(SubMat_proba[nuc]) for sub_proba in SubMat_proba[nuc]]

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
        FocalSeq = FocalSeqs[ID].seq
        AncestralSeq = AncestralSeqs[ID].seq

        NbSub = get_sub_number(FocalSeq, AncestralSeq)

        if NbSub > 1:
            DeltaObs = calculate_delta_svm(AncestralSeq, FocalSeq, KmerLen, SVM_dict)

            RandomSeqs = get_random_seqs(AncestralSeq, SubMat_proba, SubMat_proba_normed, SubMat_direction, NbSub, NbRand)
            DeltaRand = [calculate_delta_svm(AncestralSeq, RandSeq, KmerLen, SVM_dict) for RandSeq in RandomSeqs]

            # Calculate p-value
            NbHigherRand = sum(rand > DeltaObs for rand in DeltaRand)
            pvalHigh = NbHigherRand / len(DeltaRand)

            Output.write(ID + "\t" + str(DeltaObs) + "\t" + str(NbSub) + "\t" + str(pvalHigh) + "\n")

Output.close()

print("Done!")

####################################################################################################
