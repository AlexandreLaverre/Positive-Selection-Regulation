#!/usr/bin/env python
# coding=utf-8

import argparse
from Bio import SeqIO
from alive_progress import alive_bar
import multiprocessing.pool
import sys
from pathlib import Path

########################################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("sample", help="Study name and Transcription Factor: Wilson/CEBPA Schmidt/CTCF ...")
parser.add_argument("--peakType", default="NarrowPeaks", help="NarrowPeaks or BroadPeaks")
parser.add_argument("-T", "--NbThread", default=1, type=int, help="Number of threads for parallelization (default=8)")
args = parser.parse_args()
maxLen = 1000

path = Path(__file__).resolve().parent[4]
pathResults = f"{path}/results/positive_selection/{args.peakType}/{args.species}/{args.sample}"
sys.path.append(f"{path}/scripts/2_selection_tests/scripts/functions/")
import SVM


def run_svm_pos(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    svm_pos = []
    for pos in range(len(focal_seq)):
        svm_pos.append(str(SVM.calculate_svm_per_pos(focal_seq, pos, SVM_dict)))

    return seq_name, svm_pos


########################################################################################################################
# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(f"{pathResults}/Model/kmer_predicted_weight.txt")

# Get sequences
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(f"{pathResults}/sequences/filtered_focal_sequences.fa"), "fasta"))
SeqIDs = FocalSeqs.keys()

output = open(f"{pathResults}/SVM_per_base.txt", "w")
all_pos = '\t'.join([f"pos{i}" for i in range(maxLen)])
output.write(f"ID\t{all_pos}\n")

# Running and writing results
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            for svm in pool.imap_unordered(run_svm_pos, SeqIDs):
                bar()
                output.write(str(svm[0]) + '\t' + '\t'.join(svm[1]) + '\n')

output.close()
########################################################################################################################
