#!/usr/bin/env python
# coding=utf-8
import numpy as np
import os
import pandas
from Bio import SeqIO
from Bio.Seq import Seq
from alive_progress import alive_bar

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
pathSubMat = path + "/substitution_matrix/human/"
Focal_fasta = path + "positive_selection/human/sequences/filtered_focal_sequences.fa"

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

# Get Focal sequences
AllSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))
FirstSeqs = dict(list(AllSeqs.items())[:10000])

with open(path + "positive_selection/human/sequences/first_focal_sequences.fa", 'w') as First:
    SeqIO.write(FirstSeqs.values(), First, 'fasta')

Mutations = list(range(2, 10, 1)) + list(range(10, 110, 10))

for NbMut in Mutations:
    Simulated_fasta = open(path + "positive_selection/human/sequences/simulated_focal_sequences_" + str(NbMut) + ".fa", 'w')
    print("Running with", NbMut, "mutations per sequence...")
    with alive_bar(10000) as bar:
        for ID in FirstSeqs.keys():
            bar()  # print progress bar
            focal_seq = str(FirstSeqs[ID].seq)
            chromosome = ID.split(':')[0]

            sub_mat_proba = SubMats[chromosome]
            sub_mat_proba_norm = SubMats_norm[chromosome]

            # Mutations step
            for _ in range(NbMut):
                # Get substitution probabilities for each position
                pos_proba = np.empty(len(focal_seq))
                for pos in range(len(focal_seq)):
                    nuc = focal_seq[pos]                               # the probability to draw a given position is equal to the sum
                    pos_proba[pos] = sum(sub_mat_proba[nuc].values())  # of the substitution probabilities from this nucleotide

                normed_pos_proba = pos_proba / sum(pos_proba)  # normalisation of all positions to sum at 1

                # draw position
                rand_pos = int(np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, size=1))
                # draw directions
                old_nuc = focal_seq[rand_pos]
                directions = list(sub_mat_proba_norm[old_nuc].keys())
                proba = list(sub_mat_proba_norm[old_nuc].values())
                new_nuc = np.random.choice(directions, p=proba)
                rand_seq = focal_seq[:rand_pos] + new_nuc + focal_seq[rand_pos + 1:]  # change value in pos

                focal_seq = rand_seq

            FirstSeqs[ID].seq = Seq(focal_seq)

    SeqIO.write(FirstSeqs.values(), Simulated_fasta, 'fasta')
    Simulated_fasta.close()


