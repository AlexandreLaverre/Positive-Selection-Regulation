#!/usr/bin/env python
# coding=utf-8

import argparse
from Bio import SeqIO
from alive_progress import alive_bar
import multiprocessing.pool
from Positive_Selection_Tests import SVM_functions as SVM

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("sample", help="Study name: Wilson Schmidt ...")
parser.add_argument("TF", help="Transcription Factor name: CEBPA CTCF ...")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--sister", default=False, action='store_true', help="Run on sister's sequences instead of focal.")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
args = parser.parse_args()

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

focal_species = "sister" if args.sister else "focal"
pathSelection = f"{path}/positive_selection/{args.species}/{args.sample}/{args.TF}/"
Ancestral_fasta = pathSelection + "/sequences/filtered_ancestral_sequences.fa"
Focal_fasta = pathSelection + "/sequences/filtered_" + focal_species + "_sequences.fa"
ModelEstimation = pathSelection + "/Model/kmer_predicted_weight.txt"

Output_all = open(pathSelection + focal_species + "_all_possible_deltaSVM.txt", "w")
Output_obs = open(pathSelection + focal_species + "_observed_deltaSVM.txt", "w")


####################################################################################################
# Get deltaSVM for all possible substitutions.
def compute_all_delta(seq):
    nuc = ["A", "T", "C", "G"]
    deltas = {}
    for position in range(len(seq)):
        old_nuc = seq[position]
        for new_nuc in nuc:
            if new_nuc != old_nuc:
                ID = "pos" + str(position) + ":" + old_nuc + "-" + new_nuc
                test_seq = list(seq)
                test_seq[position] = new_nuc
                test_seq = "".join(test_seq)

                deltas[ID] = str(SVM.calculate_delta(seq, test_seq, SVM_dict))

    return deltas


# Return all and observed deltaSVM for a given sequence
def run_deltas(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)
    substitutions = SVM.get_sub_ids(ancestral_seq, focal_seq)

    if len(substitutions) > 1 and 20 <= len(ancestral_seq) <= 1000:
        SVM_score = SVM.calculate_svm(focal_seq)
        deltaSVM = SVM.calculate_delta(ancestral_seq, focal_seq, SVM_dict)

        # all possible substitutions
        deltas = compute_all_delta(ancestral_seq)
        obs_delta = '\t'.join([deltas[sub] for sub in substitutions])
        all_delta = '\t'.join(deltas.values())

        output_obs = f"{seq_name}\t{SVM_score}\t{deltaSVM}\t{len(substitutions)}\t{obs_delta}\n"
        output_all = f"{seq_name}\t{all_delta}\n"

        return output_obs, output_all


####################################################################################################
# Datas
# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(ModelEstimation)

# Get Ancestral sequences
AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))
if len(AncestralSeqs) == 0:
    raise ValueError("Ancestral sequence file is empty!")

# Get Focal sequences
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))
SeqIDs = FocalSeqs.keys()
if len(FocalSeqs) == 0:
    raise ValueError("Focal sequence file is empty!")

####################################################################################################
# Running and writing results
#Output_obs.write("ID\tSVM\tdeltaSVM\tmed.deltaSVM.simul\tmean.deltaSVM.simul\tNbSub\tpval.high\n")  # header

# protect the entry point
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(run_deltas, SeqIDs):
                bar()  # print progress bar
                if results is not None:
                    Output_obs.write(results[0])
                    Output_all.write(results[1])

Output_obs.close()
Output_all.close()
####################################################################################################
