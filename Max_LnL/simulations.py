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
from mpl_toolkits.mplot3d import Axes3D

# Seed numpy random generator
np.random.seed(1234)

# Define parameters
path = "../../results/MaxLikelihoodApproach/ProbaFixEstimation/"
os.makedirs(path, exist_ok=True)
num_simulations = 1000
n_bins = 50
max_mut = 10
distrib = "Gaussian"  # or Beta
plots = False
output = True

if distrib == "Gaussian":
    import MLEvol_functions_Gaussian as ML
else:
    import MLEvol_functions as ML

########################################################################################################################
dfs = []
for Nsim in range(num_simulations):
    print(f'Simulation {Nsim}')
    # SVM score distribution: affinity of all possible deltas for a sequence
    All_SVM = np.random.normal(-1, 3, 3000)
    gaussian_mutation = stats.gaussian_kde(All_SVM)  # fit a kernel distribution on our data
    hist_SVM = np.histogram(All_SVM, bins=n_bins)

    stab_weights = np.array([stats.norm.pdf(svm, loc=0, scale=1) for svm in All_SVM])
    stab_weights /= np.sum(stab_weights)

    pos_weights = np.array([stats.norm.pdf(svm, loc=5, scale=1) for svm in All_SVM])
    pos_weights /= np.sum(pos_weights)

    # Observed deltaSVM: change in affinity of substitutions for a sequence
    Nsub = np.random.randint(2, max_mut)
    print(f'Number of substitutions: {Nsub}')
    Simul_Obs_SVM = {"Random": np.random.choice(All_SVM, Nsub, replace=False),
                     "Stabilizing": np.random.choice(All_SVM, Nsub, replace=False, p=stab_weights),
                     "Positive": np.random.choice(All_SVM, Nsub, replace=False, p=pos_weights)}

    for simul in Simul_Obs_SVM.keys():
        print(f'Simulation {simul} {Nsim}')
        obs_svm = Simul_Obs_SVM[simul]

        # Null model: no param.
        LL_neutral = ML.loglikelihood(obs_svm, [], hist_SVM)
        print(f"  LnL neutral    : {LL_neutral:.2f}")

        # Stabilizing selection: alpha = beta
        model_purif = minimize(lambda theta: -ML.loglikelihood(obs_svm, theta, hist_SVM), np.array([1.0]),
                               bounds=[(0.0, np.inf)], method="Nelder-Mead")
        LL_purif = -model_purif.fun
        print(f"  LnL stabilizing: {LL_purif:.4f}, alpha: {model_purif.x[0]:.3g}, beta: {model_purif.x[0]:.3g}")

        bounds = [(0.0, np.inf), (0.0, np.inf)]
        # Positive selection
        initial_guess = np.array([1.0, 1.0])
        model_pos = minimize(lambda theta: -ML.loglikelihood(obs_svm, theta, hist_SVM), initial_guess,
                             bounds=bounds, method="Nelder-Mead")
        LL_pos = -model_pos.fun
        print(f"  LnL positive   : {LL_pos:.4f}, alpha: {model_pos.x[0]:.3g}, beta: {model_pos.x[1]:.3g}")

        # Likelihood ratio test:
        LRT_null_purif = -2 * (LL_neutral - LL_purif)
        LRT_purif_pos = -2 * (LL_purif - LL_pos)

        p_value_null_purif = chi2.sf(LRT_null_purif, 1)
        p_value_purif_pos = chi2.sf(LRT_purif_pos, 1)
        print(f"  LRT null vs purif: {LRT_null_purif:.2f}, p-value: {p_value_null_purif:.2g}")
        print(f"  LRT purif vs pos : {LRT_purif_pos:.2f}, p-value: {p_value_purif_pos:.2g}")

        alpha = 0.05
        conclusion = "Neutral model"
        if p_value_null_purif < alpha:
            conclusion = "Stabilizing model"
        if p_value_purif_pos < p_value_null_purif and p_value_purif_pos < alpha:
            conclusion = "Positive model"

        print("  Conclusion:", conclusion)

        if output:
            # Create a DataFrame for each simulation
            df = pd.DataFrame({"Nmut": [len(obs_svm)], "MeanObs": [np.mean(obs_svm)], "VarObs": [np.var(obs_svm)],
                               "Scenario": [simul], "AlphaPurif:": [model_purif.x[0]], "AlphaPos:": [model_pos.x[0]],
                               "BetaPos:": [model_pos.x[1]], "NiterPurif:": [model_purif.nit],
                               "NiterPos:": [model_pos.nit],
                               "LL_neutral": [LL_neutral], "LL_purif": [LL_purif], "LL_pos": [LL_pos],
                               "LRT_null_purif": [LRT_null_purif], "LRT_purif_pos": [LRT_purif_pos],
                               "p_value_null_purif": [p_value_null_purif], "p_value_purif_pos": [p_value_purif_pos],
                               "Conclusion": [conclusion]})

            dfs.append(df)

        if plots:
            with PdfPages(f'{path}/MLE_summary_{simul}_{Nsim}_{distrib}.pdf') as pdf:
                fig, axes = plt.subplots(2, 2, figsize=(14, 14))
                axes[1, 1].axis('off')
                axes[1, 1] = fig.add_subplot(224, projection='3d')
                ML.general_plot(All_SVM, obs_svm, gaussian_mutation, model_purif, model_pos, simul, conclusion,
                                axes[0, 0], axes[0, 1])
                ML.plot_model(obs_svm, hist_SVM, model_purif.x, axes[1, 0], model_type="Stabilizing", bounds=bounds)
                ML.plot_model(obs_svm, hist_SVM, model_pos.x, axes[1, 1], model_type="Positive", bounds=bounds)
                pdf.savefig()
                plt.close(fig)
                plt.clf()

if output:
    # Concatenate all individual DataFrames
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'{path}/MLE_summary_{num_simulations}simul_{n_bins}bins_{max_mut}MaxMut_{distrib}.csv', index=False)

########################################################################################################################
