#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from alive_progress import alive_bar
import multiprocessing.pool
import warnings
import argparse
from mpl_toolkits.mplot3d import Axes3D
import sys
sys.path.append('/Users/alaverre/Documents/Detecting_positive_selection/scripts/Positive_Selection_Tests/Max_LnL_Test/')
import MLEvol_functions as ML

warnings.filterwarnings("ignore", category=RuntimeWarning)


########################################################################################################################
def estimate_evolution(id, plots=False):
    # Retrieve the corresponding deltas for each ID
    all_svm_row = All_SVM_All_seq.loc[All_SVM_All_seq['ID'] == id, 1:].iloc[0]
    obs_svm_row = Obs_SVM_All_seq.loc[Obs_SVM_All_seq['ID'] == id, 4:].iloc[0]

    # SVM score distribution: affinity of all possible deltas for a sequence
    all_svm = all_svm_row.dropna().values.tolist()
    obs_svm = obs_svm_row.dropna().values.tolist()
    hist_svm = np.histogram(all_svm, bins=args.NbBin)

    estimations, models = ML.run_estimations(hist_svm, obs_svm, alpha=0.01)
    estimations.insert(0, "ID", [id])

    if plots:
        # fit a kernel distribution on our data
        gaussian_mutation = stats.gaussian_kde(all_svm)
        with PdfPages(f'{pathResults}/MLE_summary_{id}.pdf') as pdf:
            fig, axes = plt.subplots(2, 2, figsize=(14, 14))
            axes[1, 1].axis('off')
            axes[1, 1] = fig.add_subplot(224, projection='3d')
            ML.general_plot(all_svm, obs_svm, gaussian_mutation, models[0], models[1], id, estimations["Conclusion"][0],
                            axes[0, 0], axes[0, 1])
            ML.plot_model(obs_svm, hist_svm, models[0].x, axes[1, 0], model_type="Stabilizing", bounds=models[2])
            ML.plot_model(obs_svm, hist_svm, models[1].x, axes[1, 1], model_type="Positive", bounds=models[2])
            pdf.savefig()
            plt.close(fig)
            plt.clf()

    return estimations


########################################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species common name (e.g: human dog)")
parser.add_argument("sample", help="Study name: Wilson Schmidt")
parser.add_argument("TF", help="Transcription factor (e.g: CEBPA CTCF)")
parser.add_argument("cluster", default="local", help="cluster or local")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
parser.add_argument("--NbBin", default=100, type=int, required=False, help="Number of bins for All deltasSVM per Seq")
parser.add_argument("--Simulation", type=str, required=False,
                    help="Type of simulations (e.g: by_500_rounds_neutral by_deltas_stabilising)")
args = parser.parse_args()

maxSub = 150
maxLength = 1000

if args.cluster == "cluster":
    path = "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel/results/"
else:
    path = "/Users/alaverre/Documents/Detecting_positive_selection/"

pathData = f'{path}/results/positive_selection/{args.species}/{args.sample}/{args.TF}/deltas/'
pathResults = f'{path}/results/MaxLikelihoodApproach/{args.species}/{args.TF}/simulations/'
os.makedirs(pathResults, exist_ok=True)

if args.Simulation:
    Ancestral_deltas_file = "simulated_initial_all_possible_deltaSVM.txt"
    Focal_deltas_file = f"simulated_{args.Simulation}_observed_deltaSVM.txt"
    Output_file = f"MLE_summary_simulated_{args.Simulation}_{args.NbBin}bins.csv"
else:
    Ancestral_deltas_file = "ancestral_all_possible_deltaSVM.txt"
    Focal_deltas_file = "focal_observed_deltaSVM.txt"
    Output_file = f"MLE_summary_{args.NbBin}bins.csv"

All_SVM_All_seq = pd.read_csv(f'{pathData}/{Ancestral_deltas_file}', sep='\t', header=None, names=range(1+maxLength*3))
Obs_SVM_All_seq = pd.read_csv(f'{pathData}/{Focal_deltas_file}', sep='\t', header=None, names=range(maxSub+4))

All_SVM_All_seq.columns = ['ID'] + list(All_SVM_All_seq.columns[1:])
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
                bar()  # print progress bar
                if result is not None:
                    dfs.append(result)

    # Concatenate all individual DataFrames
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'{pathResults}/{Output_file}', index=False)

########################################################################################################################
