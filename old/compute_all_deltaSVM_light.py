#!/usr/bin/env python
# coding=utf-8

import argparse
from Bio import SeqIO
from alive_progress import alive_bar
import multiprocessing.pool
import sys
import os
from pathlib import Path

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("sample", help="Study name and Transcription Factor: Wilson/CEBPA Schmidt/CTCF ...")
parser.add_argument("--peakType", default="NarrowPeaks", help="NarrowPeaks or BroadPeaks")
parser.add_argument("-N", "--node", default="ancestral",
                    help="From which node to compute deltas: ancestral, focal or sister (default=ancestral)")
parser.add_argument("-S", "--Simulation", default=False,
                    help="Get obs delta for all the simulated regimes (default=False; either 500_rounds or deltas")
parser.add_argument("-T", "--NbThread", default=1, type=int, help="Number of threads for parallelization (default=8)")
args = parser.parse_args()

maxSub = 150
maxLen = 1000

path = Path(__file__).resolve().parent.parent
pathResults = f"{path}/results/positive_selection/{args.peakType}/{args.species}/{args.sample}"
sys.path.append(f"{path}/scripts/2_selection_tests/functions/")
import SVM

####################################################################################################
output_files = {}
os.makedirs(f"{pathResults}/deltas/", exist_ok=True)
output_files['all'] = open(f"{pathResults}/deltas/{args.node}_all_possible_deltaSVM.txt", "w")
output_files['focal'] = open(f"{pathResults}/deltas/{args.node}_to_observed_deltaSVM.txt", "w")

# Get reference sequences
ReferenceSeqs = SeqIO.to_dict(SeqIO.parse(open(f"{pathResults}/sequences/filtered_{args.node}_sequences.fa"), "fasta"))

# Get focal sequences
targets = []
targets.append('focal')
FocalSeqs = {'focal': SeqIO.to_dict(SeqIO.parse(open(f"{pathResults}/sequences/filtered_focal_{args.node}_sequences.fa"), "fasta"))}
SeqIDs = FocalSeqs['focal'].keys()

# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(f"{pathResults}/Model/kmer_predicted_weight.txt")

####################################################################################################
# Return all and observed deltaSVM for a given sequence
def run_deltas(seq_name):
    reference_seq = str(ReferenceSeqs[seq_name].seq)
    output_all, output_obs = "", {}

    if 20 <= len(reference_seq) <= maxLen:
        # compute all deltas from the reference node
        deltas = SVM.compute_all_delta(reference_seq, SVM_dict)
        all_delta = '\t'.join(deltas.values())
        output_all = f"{seq_name}\t{all_delta}\n"

        focal_seq = str(FocalSeqs["focal"][seq_name].seq)
        substitutions = SVM.get_sub_ids(reference_seq, focal_seq)

        if maxSub > len(substitutions) > 1:
            svm_score = SVM.calculate_svm(focal_seq, SVM_dict)
            delta_svm = SVM.calculate_delta(reference_seq, focal_seq, SVM_dict)
            obs_delta = '\t'.join([deltas[sub] for sub in substitutions])

            output_obs['focal'] = f"{seq_name}\t{svm_score}\t{delta_svm}\t{len(substitutions)}\t{obs_delta}\n"

    return output_all, output_obs


####################################################################################################
# Write header
all_mutations = '\t'.join([f"pos{i}:{nuc}" for i in range(0, maxLen) for nuc in ["A", "T", "C", "G"]])
output_files['all'].write(f"ID\t{all_mutations}\n")

# Running and writing results
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            for results in pool.imap_unordered(run_deltas, SeqIDs):
                bar()
                if results is not None:
                    output_files['all'].write(results[0])
                    if len(targets) > 0 and len(results[1]) > 0:
                        for evol in targets:
                            output_files[evol].write(results[1][evol])

for file in output_files.values():
    file.close()
####################################################################################################
