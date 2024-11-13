import random
from Bio.Seq import Seq
import os
import pandas
import numpy as np
from scipy import stats
import itertools as it
import MLEvol as ML

all_nuc = ["A", "T", "C", "G"]


########################################################################################################################
##################################### Deltas associated functions ######################################################
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
def calculate_svm(seq, svm_dict, kmer_len=10):
    if kmer_len != 10:
        kmer_len = len(list(svm_dict.keys())[0])
    svm = 0
    for pos in range(len(seq) - kmer_len + 1):  # sliding window of kmer length
        kmer = seq[pos:pos + kmer_len]
        svm += svm_dict[kmer]  # sum of SVM for each kmer

    return round(svm, 7)


# Calculate SVM from sliding windows
def calculate_svm_per_pos(seq, pos, svm_dict, kmer_len=10):
    if kmer_len != 10:
        kmer_len = len(list(svm_dict.keys())[0])
    svm = 0
    for window in range(kmer_len):  # sliding window of kmer length
        start = pos-window
        end = start+kmer_len
        if start < 0 or end > len(seq):  # if the window is out of the sequence
            continue
        kmer = seq[start:end]
        svm += svm_dict[kmer]  # sum of SVM for each kmer

    return round(svm, 7)


# Calculate delta SVM from sliding windows
def calculate_delta(seq_ref, seq_alt, svm_dict, kmer_len=10):
    if kmer_len != 10:
        kmer_len = len(list(svm_dict.keys())[0])

    delta_svm = 0
    for pos in range(len(seq_ref) - kmer_len + 1):  # sliding window of kmer length
        kmer_ref = seq_ref[pos:pos + kmer_len]
        kmer_alt = seq_alt[pos:pos + kmer_len]
        if kmer_ref == kmer_alt:
            continue
        delta_svm += svm_dict[kmer_alt] - svm_dict[kmer_ref]  # sum of delta between sequences for each kmer

    return round(delta_svm, 7)


# Get deltaSVM for all possible substitutions.
def compute_all_delta(seq, svm_dict):
    deltas = {}
    for pos in range(len(seq)):
        old_nuc = seq[pos]
        for new_nuc in all_nuc:
            mut_id = f"pos{pos}:{new_nuc}"
            if new_nuc != old_nuc:
                mut_seq = list(seq)
                mut_seq[pos] = new_nuc
                mut_seq = "".join(mut_seq)
                deltas[mut_id] = str(calculate_delta(seq, mut_seq, svm_dict))
            else:
                deltas[mut_id] = "NA"

    return deltas


def update_deltas(seq, svm_dict, deltas, new_sub):
    new_deltas = deltas.copy()
    # Find positions around the substitution
    start = max(new_sub - 9, 0)
    end = min(new_sub + 9, len(seq)-1)
    positions = it.chain(range(start, new_sub), range(new_sub+1, end+1))
    for pos in positions:
        old_nuc = seq[pos]
        old_delta_loc = pos * 3
        for new_nuc in all_nuc:
            if new_nuc != old_nuc:
                mut_seq = list(seq)
                mut_seq[pos] = new_nuc
                mut_seq = "".join(mut_seq)
                new_deltas[old_delta_loc] = float(calculate_delta(seq, mut_seq, svm_dict))

                old_delta_loc = old_delta_loc + 1

    return new_deltas


########################################################################################################################
################################## Substitutions associated functions ##################################################
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
            id = f"pos{loc}:{pos_alt}"
            sub_ids.append(id)
        loc += 1

    return sub_ids


# Get the Substitution Matrix and the normalised one for each chromosome
def get_sub_matrix(path_matrix):
    sub_mats = {}
    sub_mats_norm = {}
    for file in os.listdir(path_matrix):
        if file.endswith('.txt'):
            chrom = file.strip('.txt')

            chrom_Table = pandas.read_table(path_matrix + file, sep=' ')
            chrom_Table.index = ['A', 'C', 'G', 'T']  # change row values
            np.fill_diagonal(chrom_Table.values, 0)  # assign 0 to diagonal
            chrom_SubMat = chrom_Table.to_dict('index')  # get dict by matrix rows
            sub_mats[chrom] = chrom_SubMat
            # get normalized dict by matrix rows
            sub_mats_norm[chrom] = chrom_Table.div(chrom_Table.sum(axis=1), axis=0).to_dict('index')

    if len(sub_mats) == 0:
        raise ValueError("Substitution matrix not found!")

    return sub_mats, sub_mats_norm


########################################################################################################################
##################################### Mutations associated functions ###################################################
def normalised_position_probability(seq, sub_prob):
    # Get substitution probabilities for each position
    pos_proba = np.empty(len(seq))
    for pos in range(len(seq)):
        nuc = seq[pos]  # the probability to draw a given position is equal to the sum
        pos_proba[pos] = sum(sub_prob[nuc].values())  # of the substitution probabilities from this nucleotide

    normed_pos_proba = pos_proba / sum(pos_proba)  # normalisation of all positions to sum at 1

    return normed_pos_proba


def normalised_mutations_probability(seq, sub_prob, sub_prob_norm):
    mut_proba = np.empty(len(seq) * 3)
    pos_mut = 0
    for pos in range(len(seq)):
        nuc = seq[pos]
        proba_pos = sum(sub_prob[nuc].values())
        proba_dir = [sub_prob_norm[nuc][i] for i in all_nuc if i != nuc]

        mut_proba[pos_mut:pos_mut + 3] = [i * proba_pos for i in proba_dir]
        pos_mut += 3

    normed_mut_proba = mut_proba / sum(mut_proba)  # normalisation of all positions to sum at 1

    return normed_mut_proba


# Mutate a sequence N time according to positional and directional probabilities
def mutate_seq_from_proba(seq, normed_pos_proba, sub_prob_norm, n_sub):
    mut_seq = seq.copy()
    # draw positions
    rand_pos = np.random.choice(np.arange(normed_pos_proba.size), p=normed_pos_proba, replace=False, size=n_sub)
    # draw directions
    for pos in rand_pos:
        old_nuc = mut_seq[pos]
        directions = list(sub_prob_norm[old_nuc].keys())
        proba = list(sub_prob_norm[old_nuc].values())
        new_nuc = np.random.choice(directions, p=proba)
        mut_seq[pos] = new_nuc

    return ''.join(mut_seq)


def get_random_seqs(seq, sub_prob, sub_prob_norm, n_sub, n_rand=1):
    normed_pos_proba = normalised_position_probability(seq, sub_prob)
    random_seqs = [""] * n_rand
    for n in range(n_rand):
        rand_seq = mutate_seq_from_proba(list(seq), normed_pos_proba, sub_prob_norm, n_sub)
        random_seqs[n] = rand_seq

    return random_seqs if n_rand != 1 else random_seqs[0]


def mutate_from_ids(seq, ids):
    ids = list(ids) if type(ids) is not list else ids
    seq = list(seq.copy())
    for sub in ids:
        id_pos = sub.split(":")[0]
        pos = int(id_pos.strip("pos"))
        old_nuc = seq[pos]
        new_nuc = sub.split(":")[1]
        seq[pos] = new_nuc
        if new_nuc == old_nuc:
            print("Warning: old and new nucleotides are the same,", sub)

    return ''.join(seq)


def mutate_from_deltas(seq, dic_deltas, n_sub, evol="random"):
    deltas = np.array(list(dic_deltas.values()))
    if evol == "random":
        sampled_sub = np.random.choice(list(dic_deltas.keys()), size=n_sub, replace=False)
    else:
        side = np.random.choice([0.99, 0.01])
        mean = np.quantile(deltas, side) if evol == "positive" else 0
        sd = 0.5 if evol == "positive" else 0.1
        weights = np.array([stats.norm.pdf(svm, loc=mean, scale=sd) for svm in deltas])
        weights /= np.sum(weights)
        sampled_sub = np.random.choice(list(dic_deltas.keys()), size=n_sub, p=weights, replace=False)

    mutate_seq = mutate_from_ids(list(seq), sampled_sub)

    return mutate_seq


def proba_delta_mut(original_seq, sub_mat, all_deltas, params, n_bins=False):
    n_deltas = len(all_deltas)
    # If random: all substitutions probabilities = 1
    if len(params) == 0:
        proba_substitution = np.ones(n_deltas)
    else:
        if n_bins:
            # Weighted probability of substitutions for each bin
            hist_svm = np.histogram(all_deltas, bins=n_bins)
            bins_values = hist_svm[1]
            proba_delta = hist_svm[0] / np.sum(hist_svm[0])
            proba_substi_bin = ML.proba_substitution(params, proba_delta, bins_values)

            # Find back the probability for each bin of deltas
            deltas_bin = np.searchsorted(bins_values, all_deltas, side='left') - 1
            proba_substitution = [proba_substi_bin[i] if i >= 0 else proba_substi_bin[0] for i in deltas_bin]

        else:
            # Compute probability of fixation for all deltas
            #delta_bounds = [np.nanmin(all_deltas), np.nanmax(all_deltas)]
            proba_fixation = np.zeros(n_deltas)
            for d in range(n_deltas):
                s = ML.coeff_selection(all_deltas[d], params)
                proba_fixation[d] = ML.proba_fixation(s)

    # Weighted probability of mutations for each position*direction (length_seq*3)
    proba_mutation = normalised_mutations_probability(original_seq, sub_mat[0], sub_mat[1])

    # Final probability = proba_mut * proba_fix
    output_array = [mut * fix for mut, fix in zip(proba_mutation, proba_fixation)]
    sum_output = np.sum(output_array)

    return list(output_array / sum_output) if sum_output != 0 else list(output_array)

########################################################################################################################
