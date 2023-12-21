#!/usr/bin/env python
# coding=utf-8
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from multiprocessing import Pool
from alive_progress import alive_bar
import sys
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/')
import SVM_functions as SVM

path = f"/Users/alaverre/Documents/Detecting_positive_selection/results/"
species = "human"
TF = "Wilson/CEBPA"
PathSequence = f"{path}/positive_selection/{species}/{TF}/sequences/"
PathModel = f"{path}/positive_selection/{species}/{TF}/Model/kmer_predicted_weight.txt"

max_mut = 10
NbThread = 8


####################################################################################################
def get_simulated_sequences(seq_id):
    original_seq = str(initial_sequences[seq_id].seq)
    if 20 <= len(original_seq) <= 1000:
        chromosome = seq_id.split(':')[0]
        sub_mat_proba = SubMats[chromosome]
        sub_mat_proba_norm = SubMats_norm[chromosome]
        nsub = np.random.randint(2, max_mut)

        # Simulate 500 sequences
        random_seq = SVM.get_random_seqs(original_seq, sub_mat_proba, sub_mat_proba_norm, n_sub=nsub, n_rand=500)
        deltas = {}
        for seq in random_seq:
            deltas[seq] = SVM.calculate_delta(original_seq, seq, SVM_dict)

        # Stabilising selection: find the lowest change in delta
        min_seq = min(deltas, key=lambda k: abs(deltas[k]))
        stab_seq = SeqRecord(Seq(min_seq), id=seq_id, description="")

        # Positive selection: find the highest change in delta =
        max_seq = max(deltas, key=lambda k: abs(deltas[k]))
        pos_seq = SeqRecord(Seq(max_seq), id=seq_id, description="")

        # Random drift: random seq
        rand_seq = random.choice(list(deltas.keys()))
        null_seq = SeqRecord(Seq(rand_seq), id=seq_id, description="")

        return seq_id, stab_seq, pos_seq, null_seq


####################################################################################################
# Model estimation for each kmer
SVM_dict = SVM.get_svm_dict(PathModel)

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(f"{path}/substitution_matrix/{species}/")

# Get initial sequences
initial_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_focal_sequences.fa"), "fasta"))
ancestral_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_ancestral_sequences.fa"), "fasta"))

# Find 1000 sequences with more than 1 substitution (for all deltas)
seq_ids = []
for ID in initial_sequences.keys():
    nb_sub = SVM.get_sub_number(ancestral_sequences[ID], initial_sequences[ID])
    if nb_sub > 1:
        seq_ids.append(ID)

    if len(seq_ids) == 1001:
        break

Stabilised_dict, Positive_dict, Neutral_dict = {}, {}, {}

####################################################################################################
if __name__ == '__main__':
    with alive_bar(len(seq_ids)) as bar:  # progress bar
        with Pool(NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(get_simulated_sequences, seq_ids):
                bar()
                if results is not None:  # can be None if sequence length is out-limit
                    Stabilised_dict[results[0]] = results[1]
                    Positive_dict[results[0]] = results[2]
                    Neutral_dict[results[0]] = results[3]

    dictionaries = {'stabilising': Stabilised_dict, 'positive': Positive_dict, 'neutral': Neutral_dict}

    for dict_name, dic in dictionaries.items():
        with open(f"{PathSequence}/simulated_sequences_{dict_name}_evolution_500.fa", 'w') as output:
            SeqIO.write(dic.values(), output, 'fasta')

####################################################################################################
