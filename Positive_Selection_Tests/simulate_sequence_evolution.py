#!/usr/bin/env python
# coding=utf-8
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing.pool
from alive_progress import alive_bar
from Positive_Selection_Tests import SVM_functions as SVM

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species = "human"
pathOutput = f"{path}/"

pathSimulation = path + "positive_selection/human/simulation/"
PathSequence = pathSimulation + "/sequences/filtered_focal_sequences.fa"
ModelEstimation = pathSimulation + "/Model/kmer_predicted_weight.txt"
max_mut = 10
NbThread = 8


####################################################################################################
def get_simulated_sequences(seq_id):
    seq = str(initial_sequences[seq_id].seq)
    chromosome = seq_id.split(':')[0]
    sub_mat_proba = SubMats[chromosome]
    sub_mat_proba_norm = SubMats_norm[chromosome]
    nsub = np.random.randint(2, max_mut)

    if 20 <= len(seq) <= 1000:
        null_seq = pos_seq = stab_seq = None
        attempt = 0
        while null_seq is None and pos_seq is None and stab_seq is None:
            random_seq = SVM.get_random_seqs(seq, sub_mat_proba, sub_mat_proba_norm, n_sub=nsub, n_rand=1)
            delta = SVM.calculate_delta(seq, random_seq, SVM_dict)

            # Low delta = Stabilising selection
            if -0.5 < delta < 0.5 and stab_seq is not None:
                stab_seq = random_seq

            # High shift in delta = Positive selection
            elif delta > 5 and pos_seq is not None:
                pos_seq = random_seq

            # Other = Random drift
            elif null_seq is not None:
                null_seq = random_seq

            attempt += 1
            if attempt > 1000:
                print(seq_id, "is highly constrained: skipped!")
                break

        Stabilised[seq_id] = stab_seq
        Positive[seq_id] = pos_seq
        Neutral[seq_id] = null_seq


####################################################################################################
# Model estimation for each kmer
SVM_dict = SVM.get_svm_dict(ModelEstimation)

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(f"{path}/substitution_matrix/{species}/")

# Get initial sequences
initial_sequences = SeqIO.to_dict(SeqIO.parse(open(PathSequence), "fasta"))
seq_ids = list(initial_sequences.keys())

####################################################################################################
if __name__ == '__main__':
    Stabilised, Positive, Neutral = {}, {}, {}
    with alive_bar(len(seq_ids)) as bar:  # progress bar
        with multiprocessing.Pool(NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(get_simulated_sequences, seq_ids):
                bar()

    dictionaries = {'stabilising': Stabilised, 'positive': Positive, 'neutral': Neutral}

    for dict_name, dic in dictionaries.items():
        with open(f"{pathOutput}/simulated_sequences_{dict_name}_evolution.fa", 'w') as output:
            SeqIO.write(dic.values(), output, 'fasta')

####################################################################################################
