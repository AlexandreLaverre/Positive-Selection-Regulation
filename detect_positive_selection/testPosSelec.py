#!/usr/bin/env python
# coding=utf-8
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas
from alive_progress import alive_bar
import multiprocessing.pool
import random

random.seed(12)
####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog")
parser.add_argument("sample", help="Species name: human dog")
parser.add_argument("NbRand", type=int, help="Number of random substitutions permutations per sequence")
parser.add_argument("Evol", default="uniform", help="Substitution model (default = uniform)")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
args = parser.parse_args()

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

pathSelection = path + "positive_selection/" + args.species + "/" + args.sample + "/"
Ancestral_fasta = pathSelection + "sequences/filtered_ancestral_sequences.fa"
Focal_fasta = pathSelection + "sequences/filtered_focal_sequences.fa"
ModelEstimation = pathSelection + "Model/kmer_predicted_weight.txt"
pathSubMat = path + "/substitution_matrix/" + args.species + "/"
Output = open(pathSelection + "PosSelTest_deltaSVM_" + str(args.NbRand) + "permutations.txt", "w")


####################################################################################################
# Functions
# Get number of substitutions per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq_ref, seq_alt))


# Calculate delta SVM from sliding windows
def calculate_delta_svm(seq_ref, seq_alt):
    delta_svm = 0
    for pos in range(len(seq_ref) - KmerLen + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos+KmerLen]
        kmer_alt = seq_alt[pos:pos+KmerLen]
        delta_svm += SVM_dict[kmer_alt] - SVM_dict[kmer_ref]  # sum of delta between sequences for each kmer
    return round(delta_svm, 7)


# Get random sequences according to substitution matrix
def get_random_seqs(seq, sub_prob, sub_prob_norm, nb_sub):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nuc = seq[pos]                                # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nuc].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of all positions to sum at 1

    random_seqs = []
    for perm in range(args.NbRand):
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


# Test positive selection for each sequence
def test_positive_selection(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)

    # Get corresponding substitution matrix
    chromosome = seq_name.split(':')[0]

    if chromosome in SubMats.keys():
        sub_mat_proba = SubMats[chromosome] if args.Evol != 'uniform' else SubMat_uniform
        sub_mat_proba_normed = SubMats_norm[chromosome] if args.Evol != 'uniform' else SubMat_uniform

        # Number of substitutions between Ancestral and Focal sequences
        nb_sub = get_sub_number(ancestral_seq, focal_seq)
        if nb_sub > 1:
            # Get observed and random deltas
            delta_obs = calculate_delta_svm(ancestral_seq, focal_seq)
            random_seqs = get_random_seqs(ancestral_seq, sub_mat_proba, sub_mat_proba_normed, nb_sub)
            delta_rand = [calculate_delta_svm(ancestral_seq, rand_seq) for rand_seq in random_seqs]

            # Calculate p-value
            nb_higher_rand = sum(rand > delta_obs for rand in delta_rand)
            p_val_high = nb_higher_rand / len(delta_rand)
            output = (seq_name + "\t" + str(delta_obs) + "\t" + str(nb_sub) + "\t" + str(p_val_high) + "\n")

            return output


####################################################################################################
# Datas
# Binding affinity values per kmer
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
if args.Evol == 'uniform':
    SubMat_uniform = {'A': {'A': 0, 'C': 0.333, 'G': 0.333, 'T': 0.334},
                      'C': {'A': 0.333, 'C': 0, 'G': 0.333, 'T': 0.334},
                      'G': {'A': 0.333, 'C': 0.333, 'G': 0, 'T': 0.334},
                      'T': {'A': 0.333, 'C': 0.333, 'G': 0.334, 'T': 0}}
else:
    # Get substitution matrix for each chromosome
    SubMats = {}
    SubMats_norm = {}
    for file in os.listdir(pathSubMat):
        if file.endswith('.txt'):
            chrom = file.strip('.txt')

            chrom_Table = pandas.read_table(pathSubMat + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']     # change row values
            np.fill_diagonal(chrom_Table.values, 0)      # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            SubMats[chrom] = chrom_SubMat

            # Get substitution probabilities of each nucleotide to sum at 1 for drawing directions
            chrom_Table_norm = chrom_Table.div(chrom_Table.sum(axis=1), axis=0)
            chrom_SubMat_norm = chrom_Table_norm.to_dict('index')
            SubMats_norm[chrom] = chrom_SubMat_norm

# Get Ancestral sequences
AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))

# Get Focal sequences
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))

####################################################################################################
# Running and writing results
Output.write("ID\tdeltaSVM\tNbSub\tpval.high\n")  # header

SeqIDs = FocalSeqs.keys()
# protect the entry point
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:

            # Run function for each sequence in parallel
            for result in pool.imap_unordered(test_positive_selection, SeqIDs):
                bar()  # print progress bar
                if result is not None:
                    Output.write(result)

Output.close()
####################################################################################################
