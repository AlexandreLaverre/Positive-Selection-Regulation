#!/usr/bin/env python
# coding=utf-8
import numpy as np
import os
import pandas
from Bio import SeqIO
from Bio.Seq import Seq
import copy
from alive_progress import alive_bar

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
pathSubMat = path + "/substitution_matrix/human/"
pathSimulation = path + "positive_selection/human/simulation/"
Focal_fasta = pathSimulation + "/sequences/filtered_focal_sequences.fa"
ModelEstimation = pathSimulation + "/Model/kmer_predicted_weight.txt"


def generate_mutated_sequence(sequence, sub_mat, sub_mat_norm, nb_mutations):
    for mut in range(nb_mutations):
        # Calculate substitution probabilities for each position
        pos_proba = np.empty(len(sequence))
        for pos in range(len(sequence)):
            nuc = sequence[pos]
            pos_proba[pos] = sum(sub_mat[nuc].values())

        normed_pos_proba = pos_proba / sum(pos_proba)  # normalisation of all positions to sum at 1

        # draw position
        rand_pos = int(np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, size=1))

        # draw directions
        old_nuc = sequence[rand_pos]
        directions = list(sub_mat_norm[old_nuc].keys())
        proba = list(sub_mat_norm[old_nuc].values())
        new_nuc = np.random.choice(directions, p=proba)
        rand_seq = sequence[:rand_pos] + new_nuc + sequence[rand_pos + 1:]  # change value in pos

        sequence = rand_seq

    return sequence


def calculate_delta_svm(seq_ref, seq_alt):
    delta_svm = 0
    for pos in range(len(seq_ref) - KmerLen + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos+KmerLen]
        kmer_alt = seq_alt[pos:pos+KmerLen]
        delta_svm += SVM_dict[kmer_alt] - SVM_dict[kmer_ref]  # sum of delta between sequences for each kmer
    return round(delta_svm, 7)

####################################################################################################
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
        SubMats_norm[chrom] = chrom_Table.div(chrom_Table.sum(axis=1), axis=0).to_dict('index')  # get normalized dict by matrix rows

# Get Focal sequences
AllSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))
FirstSeqs = dict(list(AllSeqs.items())[:10000])
First_focal_fasta = pathSimulation + "/sequences/first_focal_sequences.fa"

if not os.path.isfile(First_focal_fasta):
    with open(First_focal_fasta, 'w') as First:
        SeqIO.write(FirstSeqs.values(), First, 'fasta')

####################################################################################################
mutations = [2, 3] #list(range(2, 10, 1)) + list(range(10, 110, 10))
# Iterate over mutation values
for nb_mut in mutations:
    output = open(pathSimulation + "/sequences/simulated_sequences_" + str(nb_mut) + "mut_selection.fa", 'w')
    out_dic = open(pathSimulation + "/sequences/simulated_deltaSVM_" + str(nb_mut) + "mut_selection.txt", 'w')

    # Create a deep copy of the first sequences
    simul_seqs = copy.deepcopy(FirstSeqs)
    deltaSVM_dic = {}
    print("Running with", nb_mut, "mutations per sequence...")

    # Iterate over sequence IDs
    with alive_bar(10000) as bar:
        for i, ID in enumerate(FirstSeqs.keys()):
            bar()  # print progress bar
            focal_seq = str(FirstSeqs[ID].seq)
            chromosome = ID.split(':')[0]
            sub_mat_proba = SubMats[chromosome]
            sub_mat_proba_norm = SubMats_norm[chromosome]
            if len(focal_seq) >= 50:
                DeltaSVM = None
                while True:
                    mutated_seq = generate_mutated_sequence(focal_seq, sub_mat_proba, sub_mat_proba_norm, nb_mut)
                    DeltaSVM = calculate_delta_svm(mutated_seq, focal_seq)

                    # First 1000 sequences with positive selection
                    if i < 1000 and DeltaSVM >= 3:
                        break
                    # Then 1000 sequences with negative selection
                    elif 1000 <= i < 2000 and DeltaSVM <= -3:
                        break
                    # Remaining sequences
                    elif i >= 2000 and -1 < DeltaSVM < 1:
                        break

                deltaSVM_dic[ID] = DeltaSVM
                simul_seqs[ID].seq = Seq(mutated_seq)

    SeqIO.write(simul_seqs.values(), output, 'fasta')
    output.close()

    for ID, delta in deltaSVM_dic.items():
        out_dic.write(str(ID) + '\t' + str(delta) + '\n')
    out_dic.close()


