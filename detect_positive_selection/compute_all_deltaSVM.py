#!/usr/bin/env python
# coding=utf-8

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from alive_progress import alive_bar
import multiprocessing.pool

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("sample", help="Study name: Wilson Schmidt ...")
parser.add_argument("TF", help="Transcription Factor name: CEBPA CTCF ...")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
args = parser.parse_args()

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

pathSelection = f"{path}/positive_selection/{args.species}/{args.sample}/{args.TF}/"
Ancestral_fasta = pathSelection + "/sequences/filtered_ancestral_sequences.fa"
Focal_fasta = pathSelection + "/sequences/filtered_focal_sequences.fa"
Output_all = open(pathSelection + "all_possible_deltaSVM.txt", "w")
Output_obs = open(pathSelection + "observed_deltaSVM.txt", "w")

ModelEstimation = pathSelection + "/Model/kmer_predicted_weight.txt"

####################################################################################################
# Functions
# Get number of substitutions per sequence
def get_sub_IDs(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")

    sub_IDs = []
    loc = 0
    for pos_ref, pos_alt in zip(seq_ref, seq_alt):
        if pos_ref != pos_alt:
            ID = "pos" + str(loc) + ":" + pos_ref + "-" + pos_alt
            sub_IDs.append(ID)
        loc += 1

    return sub_IDs


# Calculate SVM from sliding windows
def calculate_svm(seq):
    svm = 0
    for pos in range(len(seq) - KmerLen + 1):   # sliding window of kmer length
        kmer = seq[pos:pos+KmerLen]
        svm += SVM_dict[kmer]   # sum of SVM for each kmer
    return round(svm, 7)


# Calculate delta SVM from sliding windows
def calculate_delta_svm(seq_ref, seq_alt):
    delta_svm = 0
    for pos in range(len(seq_ref) - KmerLen + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos+KmerLen]
        kmer_alt = seq_alt[pos:pos+KmerLen]

        if kmer_ref != kmer_alt:
            delta_svm += SVM_dict[kmer_alt] - SVM_dict[kmer_ref]  # sum of delta between sequences for each kmer
    return round(delta_svm, 7)


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

                deltas[ID] = str(calculate_delta_svm(seq, test_seq))

    return deltas

# Compute all deltaSVM
def run_deltas(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)
    substitutions = get_sub_IDs(ancestral_seq, focal_seq)

    if len(substitutions) > 1 and 20 <= len(ancestral_seq) <= 1000:
        SVM = calculate_svm(focal_seq)
        deltaSVM = calculate_delta_svm(ancestral_seq, focal_seq)

        # all possible substitutions
        deltas = compute_all_delta(ancestral_seq)
        obs_delta = '\t'.join([deltas[sub] for sub in substitutions])
        all_delta = '\t'.join(deltas.values())

        output_obs = f"{seq_name}\t{SVM}\t{deltaSVM}\t{len(substitutions)}\t{obs_delta}\n"
        output_all = f"{seq_name}\t{all_delta}\n"

        return output_obs, output_all

####################################################################################################
# Datas
# Binding affinity values per kmer
SVM_dict = {}
with open(ModelEstimation, 'r') as model:
    for i in model.readlines():
        i = i.strip("\n")
        i = i.split("\t")
        kmer = i[0]
        KmerLen = len(kmer)
        svm_score = float(i[1])
        rev_kmer = str(Seq(kmer).reverse_complement())  # add the reverse complement kmer

        SVM_dict[kmer] = svm_score
        SVM_dict[rev_kmer] = svm_score


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
