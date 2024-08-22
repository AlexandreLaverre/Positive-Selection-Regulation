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
import pandas as pd
sys.path.append("/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/functions/")
import SVM


####################################################################################################
species = "human"
TF = "Wilson/CEBPA"
path = f"/Users/alaverre/Documents/Detecting_positive_selection/results/"
pathResult = f"{path}/positive_selection/{species}/{TF}/"

max_sub = 25
NbThread = 8
maxLength = 1000
Nbin = False
epistasis = True
BackMutation = False
prefix = f"{'epistasis' if epistasis else 'independent_SVM'}_{'with_backMut' if BackMutation else 'without_backMut'}"
prefix = prefix if Nbin else f"{prefix}_noBin"


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

        sum_output = np.sum(subs_proba)
        subs_proba = list(subs_proba / sum_output) if sum_output != 0 else list(subs_proba)

    return sub_ids, subs_proba


def get_simulated_sequences(seq_id):
    original_seq = str(initial_sequences[seq_id].seq)
    seq_length = len(original_seq)
    if 20 <= seq_length <= 1000:
        chromosome = seq_id.split(':')[0]
        sub_mat = [SubMats[chromosome], SubMats_norm[chromosome]]
        nsub = np.random.randint(2, max_sub+1)

        all_svm_row = All_SVM_All_seq.loc[All_SVM_All_seq['ID'] == seq_id, "pos0:A":f"pos{seq_length-1}:G"].iloc[0]
        all_svm_no_nan = all_svm_row.dropna()
        init_sub_id = all_svm_no_nan.index.tolist()
        init_deltas = [float(svm) for svm in all_svm_no_nan.values.tolist()]

        # Beta params
        alpha_stab = np.random.randint(2000, 3000)
        pos = np.random.randint(45, 50)
        stab_params = [alpha_stab, alpha_stab]          # alpha=beta
        pos_params = random.sample([pos, 1], 2)  # alpha != beta
        all_params = {"neutral": [], "stabilising": stab_params, "positive": pos_params}

        output_proba_sub = {}
        mut_seqs = {}
        for evol in evol_regimes:
            params = all_params[evol]
            # Compute initial probabilities of substitutions (= proba_mut * proba_delta * proba_fix)
            init_proba_subs = SVM.proba_delta_mut(original_seq, sub_mat, init_deltas, params, Nbin)
            output_proba_sub[evol] = init_proba_subs

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

                    proba_subs = SVM.proba_delta_mut(mutated_seq, sub_mat, deltas, params, Nbin)

            # Mutate seq with all independent sampled substitutions at once
            if not epistasis:
                mutated_seq = SVM.mutate_from_ids(original_seq, samp_sub_ids)

            mut_seqs[evol] = mutated_seq

            substitutions = SVM.get_sub_ids(original_seq, mutated_seq)
            if not set(substitutions) == set(samp_sub_ids):
                print(seq_id, nsub, "\n sampled:", len(samp_sub_ids), set(samp_sub_ids),
                      "\n mutated:", len(substitutions), set(substitutions))

        rand_seq = SeqRecord(Seq(mut_seqs["neutral"]), id=seq_id, description="")
        stab_seq = SeqRecord(Seq(mut_seqs["stabilising"]), id=seq_id, description="")
        pos_seq = SeqRecord(Seq(mut_seqs["positive"]), id=seq_id, description="")
        sequences = {"neutral": rand_seq, "stabilising": stab_seq, "positive": pos_seq}

        stats = [str(nsub), str(alpha_stab), str(pos_params[0]), str(pos_params[1])]

        return seq_id, stats, sequences, output_proba_sub


####################################################################################################
# Model estimation for each kmer
SVM_dict = SVM.get_svm_dict(f"{pathResult}/Model/kmer_predicted_weight.txt")

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(f"{path}/substitution_matrix/{species}/")

# Get initial sequences
initial_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{pathResult}/sequences/filtered_focal_sequences.fa"), "fasta"))

# All deltaSVM
All_SVM_All_seq = pd.read_csv(f'{pathResult}/deltas/focal_all_possible_deltaSVM.txt', sep='\t', header=0)

# Find 1000 sequences with more than 1 substitution (for all deltas)
seq_ids = list(initial_sequences.keys())[0:1000]

Stats = {}
evol_regimes = ["neutral", "stabilising", "positive"]
Sequences = {evol: {} for evol in evol_regimes}
Proba_Sub = {evol: {} for evol in evol_regimes}

####################################################################################################
if __name__ == '__main__':
    with alive_bar(len(seq_ids)) as bar:  # progress bar
        with Pool(NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(get_simulated_sequences, seq_ids):
                bar()
                if results is not None:  # can be None if sequence length is out-limit
                    ID = results[0]
                    Stats[ID] = results[1]
                    for evol in evol_regimes:
                        Sequences[evol][ID] = results[2][evol]
                        Proba_Sub[evol][ID] = results[3][evol]

    # Write sequences
    for evol, seq_dict in Sequences.items():
        with open(f"{pathResult}/sequences/simul_by_params_{prefix}_{evol}_evolution.fa", 'w') as outSeq:
            for ID, seq in seq_dict.items():
                SeqIO.write(seq, outSeq, 'fasta')

    # Write Substitution probabilities
    for evol, sub_dict in Proba_Sub.items():
        with open(f"{pathResult}/substitution_probabilities/simul_by_params_{prefix}_{evol}_evolution.txt", 'w') as outSub:
            for ID, sub in sub_dict.items():
                outSub.write(ID + "\t" + "\t".join(map(str, sub)) + "\n")

    # Write stats
    with open(f"{pathResult}/stats_simulation/simul_by_params_{prefix}.txt", 'w') as outStats:
        outStats.write("ID\tNsub\tAlphaStab\tAlphaPos\tBetaPos\n")
        for ID, values in Stats.items():
            outStats.write(ID + "\t" + "\t".join(values) + "\n")

####################################################################################################
