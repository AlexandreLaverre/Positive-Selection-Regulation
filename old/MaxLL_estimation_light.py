#!/usr/bin/env python
# coding=utf-8
import numpy as np
from scipy import stats
import pandas as pd
from alive_progress import alive_bar
import multiprocessing.pool
import warnings
import argparse
import sys
from pathlib import Path

np.random.seed(1234)
warnings.filterwarnings("ignore", category=RuntimeWarning)

########################################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species common name (e.g: human dog)")
parser.add_argument("sample", help="Study and Transcription Factor: Wilson/CEBPA Schmidt/CTCF...")
parser.add_argument("--peakType", default="NarrowPeaks", help="NarrowPeaks or BroadPeaks")
parser.add_argument("--NbBin", default=50, type=int, required=False, help="Number of bins for deltasSVM (default=50)")
parser.add_argument("-N", "--node", default="ancestral", help="From which node to compute deltas (default=ancestral)")
parser.add_argument("-S", "--Simulation", default=False, help="Name of the simulation (default=False)")
parser.add_argument("-T", "--NbThread", default=8, type=int, help="Number of threads for parallelization (default=8)")
parser.add_argument("--binType", default="quantile", type=str, help="Method to cut SVM distribution (hist or quantile)")
parser.add_argument("--threshold", default="0.01", type=float, help="Significance threshold for the test (default=0.01)")
args = parser.parse_args()

maxSub = 150
maxLength = 1000

path = Path(__file__).resolve().parent.parent
pathResults = f'{path}/results/positive_selection/{args.peakType}/{args.species}/{args.sample}/'
pathMatrix = f"{path}/results/substitution_matrix/{args.species}/"

sys.path.append(f"{path}/scripts/2_selection_tests/functions/")
import MLEvol as ML
import SVM

# Get substitution matrix for each chromosome
SubMats, SubMats_norm = SVM.get_sub_matrix(pathMatrix)


########################################################################################################################
def estimate_evolution(id, plots=False):
    chromosome = id.split(':')[0] if args.sample != "CTCF_binding/CTCF" else id.split('_')[0]
    if chromosome not in SubMats:
        print(f"Chromosome {chromosome} not found in the substitution matrix")
        return None

    sub_mat_proba = SubMats[chromosome]

    # Retrieve the corresponding deltas for each ID
    all_svm_row = All_SVM_All_seq.loc[All_SVM_All_seq['ID'] == id, "pos0:A":].iloc[0]
    obs_svm_row = Obs_SVM_All_seq.loc[Obs_SVM_All_seq['ID'] == id, 4:].iloc[0]

    # SVM score distribution: affinity of all possible deltas for a sequence
    all_svm = all_svm_row.dropna().values.tolist()
    obs_svm = obs_svm_row.dropna().values.tolist()

    # Get column names of the original sequence (deltaSVM=NA)
    seq_len = int(len(all_svm)/3)
    original_seq = all_svm_row[all_svm_row.isna()].index.tolist()[0:seq_len]

    # Get all svm ids (pos:nuc:sub)
    all_svm_id = ML.get_svm_ids(original_seq)

    if len(set(obs_svm)) < 2:
        print(f"{id} has not enough substitutions with distinct effect to perform the test")
        return None

    # Run the Maximum Likelihood Estimation
    estimations, models = ML.run_estimations(all_svm=all_svm, all_svm_id=all_svm_id, obs_svm=obs_svm,
                                             sub_mat_proba=sub_mat_proba, alpha_threshold=args.threshold,
                                             min_bin=args.NbBin, bins=args.binType)
    estimations.insert(0, "ID", [id])

    return estimations


########################################################################################################################
Ancestral_deltas_file = f"{args.node}_all_possible_deltaSVM.txt"
Focal_deltas_file = f"{args.node}_to_observed_deltaSVM.txt"
Output_file = f"Tests/MLE_summary_{args.binType}_{args.node}.csv"

All_SVM_All_seq = pd.read_csv(f'{pathResults}/deltas/{Ancestral_deltas_file}', sep='\t', header=0)
Obs_SVM_All_seq = pd.read_csv(f'{pathResults}/deltas/{Focal_deltas_file}', sep='\t', header=None, names=range(maxSub+4))
Obs_SVM_All_seq.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(Obs_SVM_All_seq.columns[4:])

########################################################################################################################
SeqIDs = set(Obs_SVM_All_seq['ID'])
dfs = []

# protect the entry point
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            # Run function for each sequence in parallel
            for result in pool.imap_unordered(estimate_evolution, SeqIDs):
                if result is not None:
                    bar()  # print progress bar
                    dfs.append(result)

    # Concatenate all individual DataFrames
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'{pathResults}/{Output_file}', index=False)

########################################################################################################################


