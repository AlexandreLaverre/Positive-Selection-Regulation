#!/usr/bin/env python
# coding=utf-8

import argparse
from Bio import SeqIO
from alive_progress import alive_bar
import multiprocessing.pool
import sys
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/')
import SVM_functions as SVM

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("sample", help="Study name: Wilson Schmidt ...")
parser.add_argument("TF", help="Transcription Factor name: CEBPA CTCF ...")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--sister", default=False, action='store_true', help="Run on sister's sequences instead of focal.")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
parser.add_argument("--Simulation", default=False, action='store_true', help="Get obs delta for all the simulated regimes")
args = parser.parse_args()

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

focal_sp = "sister" if args.sister else "focal"
pathResults = f"{path}/positive_selection/{args.species}/{args.sample}/{args.TF}/"

output_files = {}
if args.Simulation:
    output_files['all'] = open(f"{pathResults}/deltas/simulated_initial_all_possible_deltaSVM.txt", "w")
    targets = ["stabilising", "neutral", "positive"]
    for evol in targets:
        output_files[evol] = open(f"{pathResults}/deltas/simulated_500_{evol}_observed_deltaSVM.txt", "w")
else:
    targets = [focal_sp]
    output_files['all'] = open(f"{pathResults}/deltas/ancestral_all_possible_deltaSVM.txt", "w")
    output_files[focal_sp] = open(f"{pathResults}/deltas/{focal_sp}_observed_deltaSVM.txt", "w")


####################################################################################################
# Get deltaSVM for all possible substitutions.
def compute_all_delta(seq):
    nuc = ["A", "T", "C", "G"]
    deltas = {}
    for position in range(len(seq)):
        old_nuc = seq[position]
        for new_nuc in nuc:
            if new_nuc != old_nuc:
                id = "pos" + str(position) + ":" + old_nuc + "-" + new_nuc
                test_seq = list(seq)
                test_seq[position] = new_nuc
                test_seq = "".join(test_seq)

                deltas[id] = str(SVM.calculate_delta(seq, test_seq, SVM_dict))

    return deltas


# Return all and observed deltaSVM for a given sequence
def run_deltas(seq_name):
    ancestral_seq = str(AncestralSeqs[seq_name].seq)
    if 20 <= len(ancestral_seq) <= 1000:
        # all possible substitutions
        deltas = compute_all_delta(ancestral_seq)
        all_delta = '\t'.join(deltas.values())
        output_all = f"{seq_name}\t{all_delta}\n"
        output_obs = {}
        for dict_name, focal in FocalSeqs.items():
            focal_seq = str(focal[seq_name].seq)
            substitutions = SVM.get_sub_ids(ancestral_seq, focal_seq)

            if len(substitutions) > 1:
                svm_score = SVM.calculate_svm(focal_seq, SVM_dict)
                delta_svm = SVM.calculate_delta(ancestral_seq, focal_seq, SVM_dict)

                obs_delta = '\t'.join([deltas[sub] for sub in substitutions])
                output_obs[dict_name] = f"{seq_name}\t{svm_score}\t{delta_svm}\t{len(substitutions)}\t{obs_delta}\n"

        return output_all, output_obs


####################################################################################################
# Datas
# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(f"{pathResults}/Model/kmer_predicted_weight.txt")

# Get sequences
if args.Simulation:
    # Get initial sequences
    AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(f"{pathResults}/sequences/filtered_focal_sequences.fa"), "fasta"))
    # Get simulated sequences
    FocalSeqs = {}
    for evol in targets:
        FocalSeqs[evol] = SeqIO.to_dict(SeqIO.parse(
            open(f"{pathResults}/sequences/simulated_sequences_{evol}_evolution_500.fa"), "fasta"))
    SeqIDs = FocalSeqs[evol].keys()
else:
    # Get ancestral sequences
    AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(f"{pathResults}/sequences/filtered_ancestral_sequences.fa"), "fasta"))
    # Get focal sequences
    FocalSeqs = {focal_sp: SeqIO.to_dict(
        SeqIO.parse(open(f"{pathResults}/sequences/filtered_{focal_sp}_sequences.fa"), "fasta"))}
    SeqIDs = FocalSeqs[focal_sp].keys()

####################################################################################################
# Running and writing results
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            for results in pool.imap_unordered(run_deltas, SeqIDs):
                bar()
                if results is not None:
                    output_files['all'].write(results[0])
                    for evol in targets:
                        output_files[evol].write(results[1][evol])

for file in output_files.values():
    file.close()
####################################################################################################
