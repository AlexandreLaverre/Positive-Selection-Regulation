#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.stats import chi2
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from alive_progress import alive_bar
import multiprocessing.pool
import warnings
import argparse
from mpl_toolkits.mplot3d import Axes3D
import MLEvol_functions as ML
warnings.filterwarnings("ignore", category=RuntimeWarning)


########################################################################################################################
def EstimateEvolution(ID, plots=False):
    #print(ID)
    # Retrieve the corresponding deltas for each ID
    all_svm_row = All_SVM_All_seq.loc[All_SVM_All_seq['ID'] == ID, 1:].iloc[0]
    obs_svm_row = Obs_SVM_All_seq.loc[Obs_SVM_All_seq['ID'] == ID, 4:].iloc[0]

    # SVM score distribution: affinity of all possible deltas for a sequence
    All_SVM = all_svm_row.dropna().values.tolist()
    obs_svm = obs_svm_row.dropna().values.tolist()

    # fit a kernel distribution on our data
    gaussian_mutation = stats.gaussian_kde(All_SVM)
    hist_SVM = np.histogram(All_SVM, bins=args.Nbins)

    # Null model: no param.
    LL_neutral = ML.loglikelihood(obs_svm, [], hist_SVM)

    std_max = np.std((min(All_SVM), max(All_SVM)))
    std_min = 1e-1
    # Stabilizing selection: fixed mean to 0, minimize variance
    model_purif = minimize(lambda theta: -ML.loglikelihood(obs_svm, theta, hist_SVM), np.array([std_max / 2]),
                           bounds=[(std_min, std_max)], method="L-BFGS-B")
    LL_purif = -model_purif.fun

    optimal_min = hist_SVM[1][1]
    optimal_max = hist_SVM[1][-2]
    bounds = [(std_min, std_max), (optimal_min, optimal_max)]
    # Positive selection: minimize mean AND variance
    model_pos = minimize(lambda theta: -ML.loglikelihood(obs_svm, theta, hist_SVM), np.array([model_purif.x[0], 0]),
                         bounds=bounds, method="Powell")
    LL_pos = -model_pos.fun

    # Likelihood ratio test:
    LRT_null_purif = -2 * (LL_neutral - LL_purif)
    LRT_null_pos = -2 * (LL_neutral - LL_pos)
    LRT_purif_pos = -2 * (LL_purif - LL_pos)

    p_value_null_purif = chi2.sf(LRT_null_purif, 1)
    p_value_null_pos = chi2.sf(LRT_null_pos, 2)
    p_value_purif_pos = chi2.sf(LRT_purif_pos, 1)

    alpha = 0.01
    if p_value_null_purif < alpha:
        conclusion = "Stabilizing model"
        if p_value_purif_pos < alpha:
            conclusion = "Positive model"
    else:
        if p_value_null_pos < alpha:
            conclusion = "Positive model"
        else:
            conclusion = "Neutral model"

    # Create a DataFrame for each simulation
    df = pd.DataFrame({"Nmut": [len(obs_svm)], "MeanObs": [np.mean(obs_svm)], "VarObs": [np.var(obs_svm)],
                       "VarPurif:": [model_purif.x[0]], "VarPos:": [model_pos.x[0]], "MeanPos:": [model_pos.x[1]],
                       "NiterPurif:": [model_purif.nit], "NiterPos:": [model_pos.nit],
                       "LL_neutral": [LL_neutral], "LL_purif": [LL_purif], "LL_pos": [LL_pos],
                       "Conclusion": [conclusion]})

    if plots:
        with PdfPages(f'{pathResults}/MLE_summary_{ID}.pdf') as pdf:
            fig, axes = plt.subplots(2, 2, figsize=(14, 14))
            axes[1, 1].axis('off')
            axes[1, 1] = fig.add_subplot(224, projection='3d')
            ML.general_plot(All_SVM, obs_svm, gaussian_mutation, model_purif, model_pos, ID, conclusion,
                            axes[0, 0], axes[0, 1])
            ML.plot_model(obs_svm, hist_SVM, model_purif.x, axes[1, 0], model_type="Stabilizing", bounds=bounds)
            ML.plot_model(obs_svm, hist_SVM, model_pos.x, axes[1, 1], model_type="Positive", bounds=bounds)
            pdf.savefig()
            plt.close(fig)
            plt.clf()

    return df


########################################################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("species", help="Species name: human dog")
parser.add_argument("TF", help="Transcription factor : CEBPA CTCF...")
parser.add_argument("--NbThread", default=1, type=int, help="Number of threads for parallelization (default = 1)")
parser.add_argument("--Nbins", default=100, type=int, required=False, help="Number of Mutation per Seq in simulation mode")
args = parser.parse_args()

maxSub = 150
maxLength = 1000

pathData = f'../../results/positive_selection/all_deltas/{args.species}/{args.TF}/'
pathResults = f'../../results/MaxLikelihoodApproach/{args.species}/{args.TF}/'
os.makedirs(pathResults, exist_ok=True)

All_SVM_All_seq = pd.read_csv(f'{pathData}/all_possible_deltaSVM.txt', sep='\t', header=None, names=range(1+maxLength*3))
All_SVM_All_seq.columns = ['ID'] + list(All_SVM_All_seq.columns[1:])

Obs_SVM_All_seq = pd.read_csv(f'{pathData}/observed_deltaSVM.txt', sep='\t', header=None, names=range(maxSub+4))
Obs_SVM_All_seq.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(Obs_SVM_All_seq.columns[4:])

SeqIDs = set(Obs_SVM_All_seq['ID'])
dfs = []
# protect the entry point
if __name__ == '__main__':
    with alive_bar(len(SeqIDs)) as bar:  # progress bar
        with multiprocessing.Pool(args.NbThread) as pool:
            # Run function for each sequence in parallel
            for result in pool.imap_unordered(EstimateEvolution, SeqIDs):
                bar()  # print progress bar
                if result is not None:
                    dfs.append(result)

    # Concatenate all individual DataFrames
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'{pathResults}/MLE_summary.csv', index=False)

########################################################################################################################
