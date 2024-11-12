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
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/functions/')
import SVM
np.random.seed(1234)

path = f"/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
species = "human"
TF = "Schmidt12/CTCF"
PathSequence = f"{path}/positive_selection/NarrowPeaks/{species}/{TF}/sequences/"
PathModel = f"{path}/positive_selection/NarrowPeaks/{species}/{TF}/Model/kmer_predicted_weight.txt"

max_mut = 20
Simu_method = "deltas"
NbThread = 8
Nsimul = 5000


####################################################################################################
def get_simulated_sequences(seq_id, method=Simu_method):
    original_seq = str(initial_sequences[seq_id].seq)
    chromosome = seq_id.split(':')[0]
    sub_mat_proba = SubMats[chromosome]
    sub_mat_proba_norm = SubMats_norm[chromosome]
    nsub = np.random.randint(2, max_mut+1)

    if method == "500_rounds":
        # Simulate 500 sequences
        random_seq = SVM.get_random_seqs(original_seq, sub_mat_proba, sub_mat_proba_norm, n_sub=nsub, n_rand=500)
        deltas = {}
        for seq in random_seq:
            deltas[seq] = SVM.calculate_delta(original_seq, seq, SVM_dict)

        # Stabilising selection: find the lowest change in delta
        stab_id = min(deltas, key=lambda k: abs(deltas[k]))

        # Positive selection: find the highest change in delta
        pos_id = max(deltas, key=lambda k: abs(deltas[k]))

        # Random drift: find a random seq
        rand_id = random.choice(list(deltas.keys()))

    elif method == "deltas":
        deltas = SVM.compute_all_delta(original_seq, SVM_dict)
        deltas = {k: v for k, v in deltas.items() if v not in [None, "NA"]}
        for k, v in deltas.items():
            deltas[k] = float(v)

        stab_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, evol="stabilising")
        pos_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, evol="positive")
        rand_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, evol="random")

    stab_seq = SeqRecord(Seq(stab_id), id=seq_id, description="")
    pos_seq = SeqRecord(Seq(pos_id), id=seq_id, description="")
    null_seq = SeqRecord(Seq(rand_id), id=seq_id, description="")

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
    chr = ID.split(':')[0]
    if chr in SubMats.keys() and 20 <= len(initial_sequences[ID]) <= 1000 and nb_sub > 1:
        seq_ids.append(ID)

    if len(seq_ids) == Nsimul:
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
        with open(f"{PathSequence}/simulated_sequences_by_{Simu_method}_{dict_name}_evolution.fa", 'w') as output:
            SeqIO.write(dic.values(), output, 'fasta')

####################################################################################################
