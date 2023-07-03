#!/usr/bin/env python
# coding=utf-8
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import scipy.stats as stats
import pandas
from alive_progress import alive_bar
import multiprocessing.pool

np.random.seed(12)

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog")
parser.add_argument("sample", help="Study name: Wilson Schmidt")
parser.add_argument("TF", help="Transcription Factor name: CEBPA CTCF")
parser.add_argument("NbRand", type=int, help="Number of random substitutions permutations per sequence")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
parser.add_argument("--Simul", type=int, required=False, help="Number of Mutation per Seq in simulation mode")
parser.add_argument("--Selection", action='store_true', default=False, help="Add a selection step during permutations (default=False)")
parser.add_argument("--Evol", required=False, default="matrix", help="Substitution model (default = matrix)")
args = parser.parse_args()

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"

if args.Simul:
    pathSelection = f"{path}positive_selection/{args.species}/simulation/{args.TF}/"
    SimulSel_flag = "_selection" if args.Selection else ""

    pathJialin = "/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/mouse/sequences/"
    seq = f"{args.TF}_{args.sample}"

    Ancestral_fasta = f"{pathJialin}{seq}_filtered_ancestor.fa"
    Focal_fasta = f"{pathJialin}{seq}_filtered_focal.fa"
    Output = open(f"{pathSelection}PosSelTest_deltaSVM_mouse_triplets_{seq}{SimulSel_flag}.txt", "w")

    #Ancestral_fasta = pathSelection + "sequences/simulated_sequences_" + str(args.Simul) + "_mut.fa"
    #Focal_fasta = pathSelection + "sequences/first_focal_sequences.fa"
    #Output = open(f"{pathSelection}PosSelTest_deltaSVM_{str(args.Simul)}_mutations{SimulSel_flag}.txt", "w")
else:
    pathSelection = f"{path}/positive_selection/{args.species}/{args.sample}/{args.TF}/"
    Ancestral_fasta = pathSelection + "sequences/filtered_ancestral_sequences.fa"
    Focal_fasta = pathSelection + "sequences/filtered_focal_sequences.fa"
    Output = open(pathSelection + "PosSelTest_deltaSVM_" + str(args.NbRand) + "permutations.txt", "w")

ModelEstimation = pathSelection + "Model/kmer_predicted_weight.txt"
pathSubMat = path + "/substitution_matrix/" + args.species + "/"


####################################################################################################
# Functions
# Get number of substitutions per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq_ref, seq_alt))


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


def mutate_seq(seq, normed_pos_proba, sub_prob_norm, sub):
    # draw positions
    rand_pos = np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, replace=False, size=sub)
    # draw directions
    for pos in rand_pos:
        old_nuc = seq[pos]
        directions = list(sub_prob_norm[old_nuc].keys())
        proba = list(sub_prob_norm[old_nuc].values())
        new_nuc = np.random.choice(directions, p=proba)
        seq[pos] = new_nuc

    return ''.join(seq)


# Get random sequences according to substitution matrix
def get_random_seqs(seq, seqID, sub_prob, sub_prob_norm, sub, selection):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nuc = seq[pos]                                # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nuc].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of all positions to sum at 1

    random_seqs = [""] * args.NbRand
    old_svm = SVMAncestral[seqID]
    nb_seq = 0
    while nb_seq < args.NbRand:
        rand_seq = mutate_seq(list(seq), normed_pos_proba, sub_prob_norm, sub)

        if not selection:
            random_seqs[nb_seq] = rand_seq
            nb_seq += 1
        else:
            delta = calculate_delta_svm(seq, rand_seq)
            new_svm = old_svm + delta
            probability = delta_distribution.cdf(new_svm)  # calculate the AUC from this SVM
            if probability > np.random.random():
                random_seqs[nb_seq] = rand_seq
                nb_seq += 1

    return random_seqs


# Test positive selection for each sequence
def test_positive_selection(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)

    # Get corresponding substitution matrix
    chromosome = seq_name.split('_')[0]  # remember to change _ by :

    if chromosome in SubMats.keys():
        sub_mat_proba = SubMats[chromosome] if args.Evol != 'uniform' else SubMat_uniform
        sub_mat_proba_normed = SubMats_norm[chromosome] if args.Evol != 'uniform' else SubMat_uniform

        # Number of substitutions between Ancestral and Focal sequences
        nb_sub = get_sub_number(ancestral_seq, focal_seq)
        if nb_sub > 1 and len(focal_seq) > 40:
            SVM = calculate_svm(focal_seq)
            # Get observed and random deltas
            delta_obs = deltaObs[seq_name]
            random_seqs = get_random_seqs(ancestral_seq, seq_name, sub_mat_proba, sub_mat_proba_normed,
                                          nb_sub, selection=args.Selection)
            delta_rand = [calculate_delta_svm(ancestral_seq, rand_seq) for rand_seq in random_seqs]

            # Calculate p-value
            nb_higher_rand = sum(rand > delta_obs for rand in delta_rand)
            p_val_high = nb_higher_rand / len(delta_rand)

            output = f"{seq_name}\t{SVM}\t{delta_obs}\t{np.median(delta_rand)}\t{np.mean(delta_rand)}" \
                     f"\t{nb_sub}\t{p_val_high}\n"

            return output


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

# Substitution matrix
if args.Evol == 'uniform':
    SubMat_uniform = {'A': {'A': 0, 'C': 0.333, 'G': 0.333, 'T': 0.334},
                      'C': {'A': 0.333, 'C': 0, 'G': 0.334, 'T': 0.333},
                      'G': {'A': 0.333, 'C': 0.334, 'G': 0, 'T': 0.333},
                      'T': {'A': 0.334, 'C': 0.333, 'G': 0.333, 'T': 0}}
else:
    # Get substitution matrix for each chromosome
    SubMats = {}
    SubMats_norm = {}
    for file in os.listdir(pathSubMat):
        if file.endswith('.txt'):
            chrom = file.strip('.txt')

            chrom_Table = pandas.read_table(pathSubMat + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']     # change row values
            np.fill_diagonal(chrom_Table.values, 0)      # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            SubMats[chrom] = chrom_SubMat

            # Get substitution probabilities of each nucleotide to sum at 1 for drawing directions
            chrom_Table_norm = chrom_Table.div(chrom_Table.sum(axis=1), axis=0)
            chrom_SubMat_norm = chrom_Table_norm.to_dict('index')
            SubMats_norm[chrom] = chrom_SubMat_norm

    if len(SubMats) == 0:
        raise ValueError("Substitution matrix not found!")


# Get Ancestral sequences
AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))
if len(AncestralSeqs) == 0:
    raise ValueError("Ancestral sequence file is empty!")

# Get Focal sequences
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))
SeqIDs = FocalSeqs.keys()
if len(FocalSeqs) == 0:
    raise ValueError("Focal sequence file is empty!")

# Calculate observed deltaSVM
deltaObs = {}
SVMAncestral = {}
for ID in SeqIDs:
    focal_seq = str(FocalSeqs[ID].seq)
    ancestral_seq = str(AncestralSeqs[ID].seq)

    # Number of substitutions between Ancestral and Focal sequences
    substitution = get_sub_number(ancestral_seq, focal_seq)
    if substitution > 1 and len(focal_seq) > 40:
        delta = calculate_delta_svm(ancestral_seq, focal_seq)
        deltaObs[ID] = delta
        svm = calculate_svm(ancestral_seq)
        SVMAncestral[ID] = svm

# Fit a probability distribution to observed deltaSVMs
if args.Selection:
    fit_params = stats.norm.fit(list(SVMAncestral.values()))
    delta_distribution = stats.norm(*fit_params)

####################################################################################################
# Running and writing results
Output.write("ID\tSVM\tdeltaSVM\tmed.deltaSVM.simul\tmean.deltaSVM.simul\tNbSub\tpval.high\n")  # header

# protect the entry point
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            # Run function for each sequence in parallel
            for result in pool.imap_unordered(test_positive_selection, SeqIDs):
                bar()  # print progress bar
                if result is not None:
                    Output.write(result)

Output.close()
####################################################################################################
