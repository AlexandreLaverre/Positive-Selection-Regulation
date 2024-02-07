#!/usr/bin/env python
# coding=utf-8
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from alive_progress import alive_bar
import random
import sys
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/')
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/Max_LnL_Test/')
import SVM_functions as SVM
import pandas as pd

####################################################################################################
path = f"/Users/alaverre/Documents/Detecting_positive_selection/results/"
species = "human"
TF = "Wilson/CEBPA"
PathSequence = f"{path}/positive_selection/{species}/{TF}/sequences/"
PathModel = f"{path}/positive_selection/{species}/{TF}/Model/kmer_predicted_weight.txt"
pathData = f'{path}/positive_selection/{species}/{TF}/deltas/'

max_mut = 10
NbThread = 8
maxSub = 150
maxLength = 1000
NbBin = 50
epistasis = True
BackMutation = False
prefix = f"{'epistasis' if epistasis else 'independent_SVM'}_{'with_backMut' if BackMutation else 'without_backMut'}"


####################################################################################################
def remove_pos(sub_ids, subs_proba, sub_locs, back_mut):
    if len(sub_locs) > 0:
        for loc in sub_locs:
            if back_mut:
                # Only remove the actual substitution
                start = loc
                end = start + 1
            else:
                # Remove the 3 possible substitutions at this position
                start = (loc//3)*3
                end = start+3

            sub_ids = sub_ids[:start] + sub_ids[end:]
            subs_proba = subs_proba[:start] + subs_proba[end:]

        subs_proba = list(subs_proba / np.sum(subs_proba))

    return sub_ids, subs_proba


def get_simulated_sequences(seq_id):
    original_seq = str(initial_sequences[seq_id].seq)
    seq_length = len(original_seq)
    if 20 <= seq_length <= 1000:
        chromosome = seq_id.split(':')[0]
        sub_mat = [SubMats[chromosome], SubMats_norm[chromosome]]
        nsub = np.random.randint(2, max_mut+1)

        all_svm_row = All_SVM_All_seq.loc[All_SVM_All_seq['ID'] == seq_id, "pos0:A":f"pos{seq_length-1}:G"].iloc[0]
        all_svm_no_nan = all_svm_row.dropna()
        init_sub_id = all_svm_no_nan.index.tolist()
        init_deltas = [float(svm) for svm in all_svm_no_nan.values.tolist()]

        # Beta params
        alpha_stab = np.random.randint(500, 3000)
        alpha_pos = np.random.randint(30, 50)

        null_params = []         # no param
        stab_params = [alpha_stab, alpha_stab]          # alpha=beta
        pos_params = random.sample([alpha_pos, 1], 2)  # alpha != beta
        all_params = [null_params, stab_params, pos_params]

        mut_seqs = []
        for params in all_params:
            # Compute initial probabilities of substitutions (= proba_mut * proba_delta * proba_fix)
            init_proba_subs = SVM.proba_delta_mut(original_seq, sub_mat, init_deltas, params, NbBin)

            # Apply substitutions on original sequence
            samp_sub_ids = []
            samp_sub_locs = []
            tmp_samp_pos = []
            mutated_seq = original_seq
            deltas = init_deltas
            proba_subs = init_proba_subs
            for _ in range(nsub):
                # Remove the mutation(s) related to the sampled sub (if not BackMutation, remove 3, else 1)
                sub_ids, mutated_proba_subs = remove_pos(init_sub_id, proba_subs, samp_sub_locs, BackMutation)
                #print(len(sub_ids), len(mutated_proba_subs))

                # Sample a substitution
                sampled_sub = np.random.choice(sub_ids, p=mutated_proba_subs)
                sub_loc = sub_ids.index(sampled_sub)
                samp_sub_ids.append(sampled_sub)
                samp_sub_locs.append(sub_loc)
                tmp_samp_pos.append(sampled_sub.split(":")[0])

                # Considering epistasis (non-independence of substitutions effect)
                if epistasis:
                    # Mutate sequence for the sampled substitution
                    mutated_seq = SVM.mutate_from_ids(mutated_seq, sampled_sub)

                    # Recompute surrounding deltas and substitution probabilities
                    sub_pos = sub_loc // 3
                    deltas = SVM.update_deltas(mutated_seq, SVM_dict, deltas, sub_pos)

                    proba_subs = SVM.proba_delta_mut(mutated_seq, sub_mat, deltas, params, NbBin)

            # Mutate seq with all independent sampled substitutions at once
            if not epistasis:
                mutated_seq = SVM.mutate_from_ids(original_seq, samp_sub_ids)

            mut_seqs.append(mutated_seq)

            substitutions = SVM.get_sub_ids(original_seq, mutated_seq)
            if not set(substitutions) == set(samp_sub_ids):
                print(seq_id, nsub, "\n sampled:", len(samp_sub_ids), set(samp_sub_ids),
                      "\n mutated:", len(substitutions), set(substitutions))

        rand, stab, pos = mut_seqs[0], mut_seqs[1], mut_seqs[2]

        rand_seq = SeqRecord(Seq(rand), id=seq_id, description="")
        stab_seq = SeqRecord(Seq(stab), id=seq_id, description="")
        pos_seq = SeqRecord(Seq(pos), id=seq_id, description="")

        return seq_id, stab_seq, pos_seq, rand_seq


####################################################################################################
# Model estimation for each kmer
SVM_dict = SVM.get_svm_dict(PathModel)

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(f"{path}/substitution_matrix/{species}/")

# Get initial sequences
initial_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_focal_sequences.fa"), "fasta"))
ancestral_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_ancestral_sequences.fa"), "fasta"))

# All deltaSVM
All_SVM_All_seq = pd.read_csv(f'{pathData}/focal_all_possible_deltaSVM.txt', sep='\t', header=0)

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
        with open(f"{PathSequence}/simulated_sequences_by_params_{prefix}_{dict_name}_evolution.fa", 'w') as output:
            SeqIO.write(dic.values(), output, 'fasta')

####################################################################################################
