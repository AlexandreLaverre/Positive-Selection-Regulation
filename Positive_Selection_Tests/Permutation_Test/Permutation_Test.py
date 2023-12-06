#!/usr/bin/env python
# coding=utf-8
import argparse
from Bio import SeqIO
import numpy as np
import scipy.stats as stats
from alive_progress import alive_bar
import multiprocessing.pool
from Positive_Selection_Tests import SVM_functions as SVM

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
    pathSelection = f"{path}positive_selection/{args.species}/simulation_mutational_steps/"
    SimulSel_flag = "_selection" if args.Selection else ""

    pathJialin = "/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/mouse/sequences/"
    seq = f"{args.TF}_{args.sample}"

    #Ancestral_fasta = f"{pathJialin}{seq}_filtered_ancestor.fa"
    #Focal_fasta = f"{pathJialin}{seq}_filtered_focal.fa"
    #Output = open(f"{pathSelection}PosSelTest_deltaSVM_mouse_triplets_{seq}{SimulSel_flag}.txt", "w")

    Ancestral_fasta = pathSelection + "sequences/simulated_sequences_" + str(args.Simul) + "_mut.fa"
    Focal_fasta = pathSelection + "sequences/first_focal_sequences.fa"
    Output = open(f"{pathSelection}PosSelTest_deltaSVM_{str(args.Simul)}_mutations{SimulSel_flag}.txt", "w")
    Distrib_simul = open(f"{pathSelection}Distrib_{str(args.Simul)}_mutations{SimulSel_flag}.txt", "w")
else:
    pathSelection = f"{path}/positive_selection/{args.species}/{args.sample}/{args.TF}/"
    Ancestral_fasta = pathSelection + "sequences/filtered_ancestral_sequences.fa"
    Focal_fasta = pathSelection + "sequences/filtered_focal_sequences.fa"
    Output = open(pathSelection + "PosSelTest_deltaSVM_" + str(args.NbRand) + "permutations_selection_NegativeDelta.txt", "w")
    NegativeSet = f"{path}/positive_selection/{args.species}/delta_negative_set/{args.TF}/PosSelTest_deltaSVM_1000permutations.txt"

ModelEstimation = pathSelection + "Model/kmer_predicted_weight.txt"
pathSubMat = path + "/substitution_matrix/" + args.species + "/"


####################################################################################################
# Functions
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
def get_random_seqs(seq, sub_prob, sub_prob_norm, sub, selection):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nuc = seq[pos]                                # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nuc].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)   # normalisation of all positions to sum at 1

    random_seqs = [""] * args.NbRand
    nb_seq = 0
    while nb_seq < args.NbRand:
        rand_seq = mutate_seq(list(seq), normed_pos_proba, sub_prob_norm, sub)

        if not selection:
            random_seqs[nb_seq] = rand_seq
            nb_seq += 1
            delta = SVM.calculate_delta(seq, rand_seq, SVM_dict, KmerLen)
            Distrib_simul.write(str(delta) + "\n")
        else:
            delta = SVM.calculate_delta(seq, rand_seq, SVM_dict, KmerLen)
            probability = 2*delta_distribution.cdf(delta)  # cumulative distribution *2
            if probability > 1:
                probability = 2-probability  # to center the maximum probability at mid distribution

            if probability > np.random.random():
                random_seqs[nb_seq] = rand_seq
                nb_seq += 1

    return random_seqs


# Test positive selection for each sequence
def test_positive_selection(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)

    # Get corresponding substitution matrix
    chromosome = seq_name.split(':')[0]  # remember to change "_" by ":" if not Jialin results.
    if chromosome in SubMats.keys():
        sub_mat_proba = SubMats[chromosome] if args.Evol != 'uniform' else SubMat_uniform
        sub_mat_proba_normed = SubMats_norm[chromosome] if args.Evol != 'uniform' else SubMat_uniform

        # Number of substitutions between Ancestral and Focal sequences
        nb_sub = SVM.get_sub_number(ancestral_seq, focal_seq)

        if nb_sub > 1 and len(focal_seq) > 40:
            focalSVM = SVM.calculate_svm(focal_seq, SVM_dict, KmerLen)
            # Get observed and random deltas
            delta_obs = SVM.calculate_delta(ancestral_seq, focal_seq, SVM_dict, KmerLen)
            random_seqs = get_random_seqs(ancestral_seq, sub_mat_proba, sub_mat_proba_normed,
                                          nb_sub, selection=args.Selection)
            delta_rand = [SVM.calculate_delta(ancestral_seq, rand_seq, SVM_dict, KmerLen) for rand_seq in random_seqs]

            # Calculate p-value
            nb_higher_rand = sum(rand > delta_obs for rand in delta_rand)
            p_val_high = nb_higher_rand / len(delta_rand)

            output = f"{seq_name}\t{focalSVM}\t{delta_obs}\t{np.median(delta_rand)}\t{np.mean(delta_rand)}" \
                     f"\t{nb_sub}\t{p_val_high}\n"

            return output


####################################################################################################
# Datas
# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(ModelEstimation)
KmerLen = SVM_dict[0]

# Get substitution matrix for each chromosome
if args.Evol == 'uniform':
    SubMat_uniform = {'A': {'A': 0, 'C': 0.333, 'G': 0.333, 'T': 0.334},
                      'C': {'A': 0.333, 'C': 0, 'G': 0.334, 'T': 0.333},
                      'G': {'A': 0.333, 'C': 0.334, 'G': 0, 'T': 0.333},
                      'T': {'A': 0.334, 'C': 0.333, 'G': 0.333, 'T': 0}}
else:
    SubMats, SubMats_norm = SVM.get_sub_matrix(pathSubMat)


# Get Ancestral sequences
AncestralSeqs = SeqIO.to_dict(SeqIO.parse(open(Ancestral_fasta), "fasta"))
if len(AncestralSeqs) == 0:
    raise ValueError("Ancestral sequence file is empty!")

# Get Focal sequences
FocalSeqs = SeqIO.to_dict(SeqIO.parse(open(Focal_fasta), "fasta"))
SeqIDs = FocalSeqs.keys()
if len(FocalSeqs) == 0:
    raise ValueError("Focal sequence file is empty!")

# Fit a probability distribution to Negative deltaSVMs
if args.Selection:
    DeltaNegative = []
    with open(NegativeSet, 'r') as negative:
        for i in negative.readlines()[1:]:
            i = i.strip("\n")
            i = i.split("\t")
            DeltaNegative.append(float(i[2]))

    fit_params = stats.norm.fit(DeltaNegative)
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
