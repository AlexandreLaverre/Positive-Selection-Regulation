#!/usr/bin/env python
# coding=utf-8
import argparse
from Bio import SeqIO
import numpy as np
from alive_progress import alive_bar
import multiprocessing.pool
import sys

np.random.seed(12)

####################################################################################################
# Variables and paths
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog")
parser.add_argument("sample", help="Study name: Wilson Schmidt...")
parser.add_argument("TF", help="Transcription Factor: CEBPA CTCF ...")
parser.add_argument("--peakType", default="NarrowPeaks", help="NarrowPeaks or BroadPeaks")
parser.add_argument("-N", "--node", default="ancestral", help="From which node to compute deltas (default=ancestral)")
parser.add_argument("--NbRand", default=10000, type=int, help="Number of random substitutions permutations per sequence (default=10k")
parser.add_argument("--Evol", required=False, default="matrix", help="Substitution model (default = matrix)")
parser.add_argument("-S", "--Simulation", required=False, help="Type of simulation (i.e: 500_rounds_stabilising or deltas_neutral)")
parser.add_argument("--cluster", action='store_true', help="Needed if run on cluster")
parser.add_argument("-T", "--NbThread", required=False, default=8, type=int, help="Number of threads for parallelization (default=8)")
args = parser.parse_args()

if args.cluster:
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/"
    pathResults = f"{path}/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/"
    pathResults = f"{path}/cluster/results/"


sys.path.append(f"{path}/scripts/Positive_Selection_Tests/functions/")
import SVM

if args.Simulation:
    pathSelection = f"{pathResults}/positive_selection/{args.peakType}/{args.species}/{args.sample}/{args.TF}/"
    Focal_fasta = f"{pathSelection}/sequences/simulated_sequences_by_{args.Simulation}_evolution.fa"
    Ancestral_fasta = f"{pathSelection}/sequences/filtered_focal_sequences.fa"
    Output = open(f"{pathSelection}/Tests/PosSelTest_deltaSVM_{args.NbRand}permutations_simulation_{args.Simulation}.txt", "w")

    # pathJialin = "/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/mouse/sequences/"
    # Ancestral_fasta = f"{pathJialin}{seq}_filtered_ancestor.fa"
    # Focal_fasta = f"{pathJialin}{seq}_filtered_focal.fa"
    # Output = open(f"{pathSelection}PosSelTest_deltaSVM_mouse_triplets_{seq}{SimulSel_flag}.txt", "w")
    # Ancestral_fasta = pathSelection + "sequences/simulated_sequences_" + str(args.Simul) + "_mut.fa"

    # pathSelection = f"{path}positive_selection/{args.species}/simulation_mutational_steps/"
    # seq = f"{args.TF}_{args.sample}"
    # Focal_fasta = pathSelection + "sequences/first_focal_sequences.fa"
    # Output = open(f"{pathSelection}PosSelTest_deltaSVM_{str(args.Simul)}_mutations.txt", "w")
    # Distrib_simul = open(f"{pathSelection}Distrib_{str(args.Simul)}_mutations.txt", "w")
else:
    pathSelection = f"{pathResults}/positive_selection/{args.peakType}/{args.species}/{args.sample}/{args.TF}/"
    Ancestral_fasta = f"{pathSelection}/sequences/filtered_{args.node}_sequences.fa"
    Focal_fasta = f"{pathSelection}/sequences/filtered_focal_{args.node}_sequences.fa"
    Output = open(f"{pathSelection}/Tests/PosSelTest_deltaSVM_{str(args.NbRand)}permutations_two_tailed_{args.node}.txt", "w")
    #NegativeSet = f"{path}/results/positive_selection/{args.peakType}/{args.species}/delta_negative_set/{args.TF}/PosSelTest_deltaSVM_1000permutations.txt"

ModelEstimation = f"{pathSelection}/Model/kmer_predicted_weight.txt"
pathSubMat = f"{pathResults}/substitution_matrix/{args.species}/"


####################################################################################################
# Get random sequences according to substitution matrix
def get_random_seqs(seq, sub_prob, sub_prob_norm, sub):
    normed_pos_proba = SVM.normalised_position_probability(seq, sub_prob)

    random_seqs = [""] * args.NbRand
    nb_seq = 0
    seq = list(seq)
    while nb_seq < args.NbRand:
        rand_seq = SVM.mutate_seq_from_proba(seq, normed_pos_proba, sub_prob_norm, sub)
        random_seqs[nb_seq] = rand_seq
        nb_seq += 1

        #delta = SVM.calculate_delta(seq, rand_seq, SVM_dict)
        #Distrib_simul.write(str(delta) + "\n")

    return random_seqs


# Test positive selection for each sequence
def test_positive_selection(seq_name):
    focal_seq = str(FocalSeqs[seq_name].seq)
    ancestral_seq = str(AncestralSeqs[seq_name].seq)

    # Get corresponding substitution matrix
    chromosome = seq_name.split(':')[0] if args.sample != "CTCF_binding" else seq_name.split('_')[0]
    if chromosome in SubMats.keys():
        sub_mat_proba = SubMats[chromosome] if args.Evol != 'uniform' else SubMat_uniform
        sub_mat_proba_normed = SubMats_norm[chromosome] if args.Evol != 'uniform' else SubMat_uniform

        # Number of substitutions between Ancestral and Focal sequences
        nb_sub = SVM.get_sub_number(ancestral_seq, focal_seq)
        if nb_sub > 1 and len(focal_seq) > 40:
            focalSVM = SVM.calculate_svm(focal_seq, SVM_dict)
            # Get observed and random deltas
            delta_obs = SVM.calculate_delta(ancestral_seq, focal_seq, SVM_dict)
            random_seqs = get_random_seqs(ancestral_seq, sub_mat_proba, sub_mat_proba_normed, nb_sub)
            delta_rand = [SVM.calculate_delta(ancestral_seq, rand_seq, SVM_dict) for rand_seq in random_seqs]

            # Calculate p-value
            nb_higher_rand = sum(rand > delta_obs for rand in delta_rand)
            nb_lower_rand = sum(rand < delta_obs for rand in delta_rand)
            p_val_high = nb_higher_rand / len(delta_rand)
            p_val_two_tailed = (2*min(nb_higher_rand, nb_lower_rand)) / len(delta_rand)

            output = f"{seq_name}\t{focalSVM}\t{delta_obs}\t{np.median(delta_rand)}\t{np.mean(delta_rand)}" \
                     f"\t{nb_sub}\t{p_val_high}\t{p_val_two_tailed}\n"

            return output


####################################################################################################
# Datas
# Binding affinity values per kmer
SVM_dict = SVM.get_svm_dict(ModelEstimation)

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

####################################################################################################
# Running and writing results
Output.write("ID\tSVM\tdeltaSVM\tmed.expected.deltaSVM\tmean.expected.deltaSVM\tNbSub\tpval.high\tpval.two.tailed\n")

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
