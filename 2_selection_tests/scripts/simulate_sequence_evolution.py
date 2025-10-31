#!/usr/bin/env python
# coding=utf-8
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from multiprocessing import Pool
from alive_progress import alive_bar
import sys
from pathlib import Path

import matplotlib.pyplot as plt
np.random.seed(1234)

parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog ...")
parser.add_argument("TF", help="Study name and Transcription Factor: Wilson/CEBPA Schmidt/CTCF ...")
parser.add_argument("Method", help="How to simulate: beta, deltas, 500_rounds")
parser.add_argument("-N", "--Nsimul", default=1000, type=int, help="Number of sequences to simulate (default=1000)")
parser.add_argument("-M", "--MaxMut", default=20, type=int, help="Number of maximum mutation (default=20)")
parser.add_argument("-T", "--NbThread", default=8, type=int, help="Number of threads for parallelization (default=8)")
parser.add_argument("--peakType", default="NarrowPeaks", help="NarrowPeaks or BroadPeaks")
parser.add_argument("--quantile", default=0.01, type=float, help="Quantile for positive selection (default=0.01)")
parser.add_argument("--PropNeutral", default=0.0, type=float, help="Proportion of Neutral Sub (default=0.0)")
parser.add_argument("--ScaleNegative", default=0.0, type=float, help="Add Negative Selection correlated with deltaSVM in Random Drift")
parser.add_argument("--SymetricNegative", default=True, help="Define if Negative Selection is symetric or not (default=True)")
parser.add_argument("--MinStab", default=25, type=float, help="Value of Alpha for Stabilising Selection (default=25)")
parser.add_argument("--MinPos", default=25, type=int, help="Value of Alpha for Directional Selection (default=25)")
parser.add_argument("--addExtreme", default=0, type=int, help="Nb of extreme mutation to add in Random Drift (default=0)")
parser.add_argument("--runOnly", default="all", type=str, help="Write only one type of simulation: stabilising, positive, neutral or all (default=all")
args = parser.parse_args()


path = Path(__file__).resolve().parents[3]
PathSequence = f"{path}/results/positive_selection/{args.peakType}/{args.species}/{args.TF}/sequences/"
PathModel = f"{path}/results/positive_selection/{args.peakType}/{args.species}/{args.TF}/Model/kmer_predicted_weight.txt"
pathMatrix = f"{path}/results/substitution_matrix/{args.species}/"
sys.path.append(f"{path}/scripts/2_selection_tests/scripts/functions/")
import SVM
import MLEvol

suffix = f"{args.Method}_{args.PropNeutral}Null_{args.MinStab}Stab_{args.MinPos}Pos_{args.ScaleNegative}Purif_addExtreme{args.addExtreme}"


####################################################################################################
def add_background_neutral(vec_proba, mut_rates, prop_neutral=0.1):
    assert np.isclose(np.sum(vec_proba), 1.0)
    assert 0 <= prop_neutral <= 1
    vec_proba *= (1 - prop_neutral)
    vec_proba += prop_neutral * mut_rates / np.sum(mut_rates)
    assert np.isclose(np.sum(vec_proba), 1.0)
    return vec_proba


def add_negative_sel(vec_pheno, sub_rates, scale=0.1, symetric=True):
    assert np.isclose(np.sum(sub_rates), 1.0)
    pheno_neg = vec_pheno < 0.5
    sub_rates[pheno_neg] *= np.exp(scale * (vec_pheno[pheno_neg]-0.5))
    if symetric:
        pheno_pos = vec_pheno > 0.5
        sub_rates[pheno_pos] *= np.exp(scale * -(vec_pheno[pheno_pos]-0.5))

    sub_rates = sub_rates / sub_rates.sum()
    assert np.isclose(np.sum(sub_rates), 1.0)
    return sub_rates


def get_simulated_sequences(seq_id, method=args.Method):
    original_seq = str(initial_sequences[seq_id].seq)
    chromosome = seq_id.split(':')[0]
    sub_mat_proba = SubMats[chromosome]
    sub_mat_proba_norm = SubMats_norm[chromosome]
    #nsub = np.random.randint(2, args.MaxMut+1)
    #nsub = np.random.choice(np.arange(2, args.MaxMut + 1, 2))
    nsub = np.random.choice([28, 29, 30])

    if method == "500_rounds":
        # Simulate 500 sequences
        random_seq = SVM.get_random_seqs(original_seq, sub_mat_proba, sub_mat_proba_norm, n_sub=nsub, n_rand=500)
        deltas = {}
        for seq in random_seq:
            deltas[seq] = SVM.calculate_delta(original_seq, seq, SVM_dict)

        # Stabilising selection: find the lowest change in delta
        stab_id = min(deltas, key=lambda k: abs(deltas[k]))

        # Positive selection: find the highest change in delta
        pos_id = max(deltas, key=lambda k: abs(deltas[k]))

        # Random drift: find a random seq
        rand_id = random.choice(list(deltas.keys()))

    elif method == "deltas":
        deltas = SVM.compute_all_delta(original_seq, SVM_dict)
        deltas = {k: v for k, v in deltas.items() if v not in [None, "NA"]}
        for k, v in deltas.items():
            deltas[k] = float(v)

        stab_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, evol="stabilising")
        pos_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, quantile=args.quantile, evol="positive")
        rand_id = SVM.mutate_from_deltas(original_seq, deltas, nsub, evol="random")

    elif method == "beta":
        deltas = SVM.compute_all_delta(original_seq, SVM_dict)
        deltas = {k: float(v) for k, v in deltas.items() if v not in [None, "NA"]}

        # Sort deltas by value
        deltas = {k: v for k, v in sorted(deltas.items(), key=lambda item: item[1])}
        deltas_pos = {k: v for k, v in deltas.items() if v >= 0}
        deltas_neg = {k: v for k, v in deltas.items() if v < 0}

        # Transform deltas to phenotype: min=0, max=1, midpoint of the distribution = 0.5 (for delta = 0)
        pheno_pos = {k: 0.5 + 0.5 * (i + 1) / (len(deltas_pos) + 1) for i, k in enumerate(deltas_pos.keys())}
        pheno = {k: 0. + 0.5 * (i + 1) / (len(deltas_neg) + 1) for i, k in enumerate(deltas_neg.keys())}
        pheno.update(pheno_pos)
        assert len(pheno) == len(deltas)
        assert min(pheno.values()) > 0.
        assert max(pheno.values()) < 1.
        mut_ids = list(deltas.keys())
        for i in range(1, len(pheno)):
            assert pheno[mut_ids[i]] > pheno[mut_ids[i - 1]]
        assert list(pheno.keys()) == list(deltas.keys())

        pos_nuc_list = [(int(i.split(":")[0].replace("pos", "")), i.split(":")[1]) for i in mut_ids]
        mut_rates = np.array([sub_mat_proba[original_seq[pos]][nuc] for pos, nuc in pos_nuc_list])
        assert not np.isnan(mut_rates).any()

        # Phenotype array
        pheno_array = np.array(list(pheno.values()))
        rand_sub_rates = MLEvol.proba_substitution(np.array([]), mut_rates, pheno_array)
        if args.ScaleNegative > 0:
            # Add negative selection
            rand_sub_rates = add_negative_sel(pheno_array, rand_sub_rates,
                                              scale=args.ScaleNegative, symetric=args.SymetricNegative)

        p_stab = args.MinStab #random.randint(args.MinStab, 10)
        stab_sub_rates = MLEvol.proba_substitution(np.array([p_stab]), mut_rates, pheno_array)
        stab_sub_rates = add_background_neutral(stab_sub_rates, mut_rates, prop_neutral=args.PropNeutral)

        p_pos = args.MinPos #random.randint(args.MinPos, 25)
        p = np.array(random.sample([p_pos, 1], k=2))
        pos_sub_rates = MLEvol.proba_substitution(p, mut_rates, pheno_array)
        pos_sub_rates = add_background_neutral(pos_sub_rates, mut_rates, prop_neutral=args.PropNeutral)

        # When selection is too strong, possible mutations can be less than 2
        nb_possible_mut = len(pos_sub_rates[pos_sub_rates > 0])
        if nb_possible_mut < 2:
            print(f"Sequence {seq_id} has not enough possible mutations")
            exit(1)
        if nb_possible_mut < nsub:
            nsub = nb_possible_mut

        # Mutate sequences
        stab_id = SVM.mutate_from_sub_rates(original_seq, mut_ids, stab_sub_rates, nsub)
        pos_id = SVM.mutate_from_sub_rates(original_seq, mut_ids, pos_sub_rates, nsub)
        rand_id = SVM.mutate_from_sub_rates(original_seq, mut_ids, rand_sub_rates, nsub)
        if args.addExtreme > 0:
            mut = random.choice([mut_ids[0], mut_ids[-1]])
            rand_id = SVM.mutate_from_ids(rand_id, [mut])
    else:
        print("Method not recognized")
        exit(1)

    stab_seq = SeqRecord(Seq(stab_id), id=seq_id, description=f"alpha={p_stab}")
    pos_seq = SeqRecord(Seq(pos_id), id=seq_id, description=f"alpha={p_pos}")
    null_seq = SeqRecord(Seq(rand_id), id=seq_id, description="")

    return seq_id, stab_seq, pos_seq, null_seq


####################################################################################################
# Model estimation for each kmer
SVM_dict = SVM.get_svm_dict(PathModel)

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(pathMatrix)

# Get initial sequences
initial_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_focal_sequences.fa"), "fasta"))
ancestral_sequences = SeqIO.to_dict(SeqIO.parse(open(f"{PathSequence}/filtered_ancestral_sequences.fa"), "fasta"))

# Randomly select N sequences with already computed deltas
eligible_ids = []
for ID in initial_sequences.keys():
    nb_sub = SVM.get_sub_number(ancestral_sequences[ID], initial_sequences[ID])
    chr = ID.split(':')[0]
    if chr in SubMats.keys() and 20 <= len(initial_sequences[ID]) <= 1000 and nb_sub > 1:
        eligible_ids.append(ID)

seq_ids = np.random.choice(eligible_ids, args.Nsimul, replace=False)

Stabilised_dict, Positive_dict, Neutral_dict = {}, {}, {}

####################################################################################################
if __name__ == '__main__':
    with alive_bar(len(seq_ids)) as bar:  # progress bar
        with Pool(args.NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(get_simulated_sequences, seq_ids):
                bar()
                if results is not None:  # can be None if sequence length is out-limit
                    Stabilised_dict[results[0]] = results[1]
                    Positive_dict[results[0]] = results[2]
                    Neutral_dict[results[0]] = results[3]

    dictionaries = {'stabilising': Stabilised_dict, 'positive': Positive_dict, 'neutral': Neutral_dict}

    #with open(f"{PathSequence}/simulated_sequences_by_{args.Method}_{args.quantile}_positive_evolution.fa", 'w') as output:
    #    SeqIO.write(dictionaries["positive"].values(), output, 'fasta')

    for dict_name, dic in dictionaries.items():
        if args.runOnly == "all" or dict_name in args.runOnly:
            with open(f"{PathSequence}/simulated_sequences_by_{suffix}_{dict_name}_evolution.fa", 'w') as output:
                SeqIO.write(dic.values(), output, 'fasta')

####################################################################################################
