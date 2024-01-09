import random
from Bio.Seq import Seq
import os
import pandas
import numpy as np
from scipy import stats


# Get number of substitutions per sequence
def get_sub_number(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq_ref, seq_alt))


# Get ID of substitutions in a sequence
def get_sub_ids(seq_ref, seq_alt):
    if len(seq_ref) != len(seq_alt):
        raise ValueError("Focal and ancestral sequences don't have the same length!")
    sub_ids = []
    loc = 0
    for pos_ref, pos_alt in zip(seq_ref, seq_alt):
        if pos_ref != pos_alt:
            id = "pos" + str(loc) + ":" + pos_ref + "-" + pos_alt
            sub_ids.append(id)
        loc += 1

    return sub_ids


# Get SVM score for each kmer of a given model
def get_svm_dict(path_model):
    svm_dict = {}
    with open(path_model, 'r') as model:
        for i in model.readlines():
            i = i.strip("\n")
            i = i.split("\t")
            kmer = i[0]
            svm_score = float(i[1])
            rev_kmer = str(Seq(kmer).reverse_complement())  # add the reverse complement kmer

            svm_dict[kmer] = svm_score
            svm_dict[rev_kmer] = svm_score

    return svm_dict


# Calculate SVM from sliding windows
def calculate_svm(seq, svm_dict):
    kmer_len = len(list(svm_dict.keys())[0])
    svm = 0
    for pos in range(len(seq) - kmer_len + 1):   # sliding window of kmer length
        kmer = seq[pos:pos+kmer_len]
        svm += svm_dict[kmer]   # sum of SVM for each kmer

    return round(svm, 7)


# Calculate delta SVM from sliding windows
def calculate_delta(seq_ref, seq_alt, svm_dict, kmer_len=10):
    if kmer_len != 10:
        kmer_len = len(list(svm_dict.keys())[0])

    delta_svm = 0
    for pos in range(len(seq_ref) - kmer_len + 1):   # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos + kmer_len]
        kmer_alt = seq_alt[pos:pos + kmer_len]
        if kmer_ref == kmer_alt:
            continue
        delta_svm += svm_dict[kmer_alt] - svm_dict[kmer_ref]  # sum of delta between sequences for each kmer

    return round(delta_svm, 7)


# Get the Substitution Matrix and the normalised one for each chromosome
def get_sub_matrix(path_matrix):
    submats = {}
    submats_norm = {}
    for file in os.listdir(path_matrix):
        if file.endswith('.txt'):
            chrom = file.strip('.txt')

            chrom_Table = pandas.read_table(path_matrix + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']  # change row values
            np.fill_diagonal(chrom_Table.values, 0)  # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            submats[chrom] = chrom_SubMat
            # get normalized dict by matrix rows
            submats_norm[chrom] = chrom_Table.div(chrom_Table.sum(axis=1), axis=0).to_dict('index')

    if len(submats) == 0:
        raise ValueError("Substitution matrix not found!")

    return submats, submats_norm


def normalised_position_probability(seq, sub_prob):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nuc = seq[pos]  # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nuc].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)  # normalisation of all positions to sum at 1

    return normed_pos_proba


# Mutate a sequence N time according to positional and directional probabilities
def mutate_seq(seq, normed_pos_proba, sub_prob_norm, n_sub):
    # draw positions
    rand_pos = np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, replace=False, size=n_sub)
    # draw directions
    for pos in rand_pos:
        old_nuc = seq[pos]
        directions = list(sub_prob_norm[old_nuc].keys())
        proba = list(sub_prob_norm[old_nuc].values())
        new_nuc = np.random.choice(directions, p=proba)
        seq[pos] = new_nuc

    return ''.join(seq)


def get_random_seqs(seq, sub_prob, sub_prob_norm, n_sub, n_rand=1):
    normed_pos_proba = normalised_position_probability(seq, sub_prob)
    random_seqs = [""] * n_rand
    for n in range(n_rand):
        rand_seq = mutate_seq(list(seq), normed_pos_proba, sub_prob_norm, n_sub)
        random_seqs[n] = rand_seq

    return random_seqs if n_rand != 1 else random_seqs[0]


# Get deltaSVM for all possible substitutions.
def compute_all_delta(seq, svm_dict):
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

                deltas[id] = str(calculate_delta(seq, test_seq, svm_dict))

    return deltas


def mutate_from_deltas(seq, dic_deltas, n_sub, evol="random"):
    deltas = np.array(list(dic_deltas.values()))
    if evol == "random":
        sampled_sub = np.random.choice(list(dic_deltas.keys()), size=n_sub, replace=False)
    else:
        mean = np.quantile(deltas, 0.99) if evol == "positive" else 0
        weights = np.array([stats.norm.pdf(svm, loc=mean, scale=0.1) for svm in deltas])
        weights /= np.sum(weights)
        sampled_sub = np.random.choice(list(dic_deltas.keys()), size=n_sub, p=weights, replace=False)

    seq = list(seq)
    for sub in sampled_sub:
        id_pos = sub.split(":")[0]
        pos = int(id_pos.strip("pos"))
        new_nuc = sub.split("-")[1]
        seq[pos] = new_nuc

    return ''.join(seq)



