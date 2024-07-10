#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

# Seed numpy random generator
np.random.seed(1234)

# Define parameters
path = "../../results/MaxLikelihoodApproach/ProbaFixEstimation/"
os.makedirs(path, exist_ok=True)
num_simulations = 3
n_bins = 100
max_mut = 10
distrib = "Gaussian"  # or Beta
plots = True
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
    #hist_SVM = np.histogram(All_SVM, bins=n_bins)

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

        estimations, models = ML.run_estimations(All_SVM, obs_svm, alpha_threshold=0.01, verbose=True)
        estimations["Scenario"] = [simul]

        if output:
            dfs.append(estimations)

        if plots:
            gaussian_mutation = stats.gaussian_kde(All_SVM)  # fit a kernel distribution on our data
            with PdfPages(f'{path}/MLE_summary_{simul}_{Nsim}_{distrib}.pdf') as pdf:
                fig, axes = plt.subplots(2, 2, figsize=(14, 14))
                axes[1, 1].axis('off')
                axes[1, 1] = fig.add_subplot(224, projection='3d')
                ML.general_plot(All_SVM, obs_svm, gaussian_mutation, models[0], models[1], simul,
                                estimations["Conclusion"][0], axes[0, 0], axes[0, 1])
                ML.plot_model(obs_svm, All_SVM, models[0].x, axes[1, 0], model_type="Stabilizing", bounds=models[2])
                ML.plot_model(obs_svm, All_SVM, models[1].x, axes[1, 1], model_type="Positive", bounds=models[2])
                pdf.savefig()
                plt.close(fig)
                plt.clf()

if output:
    # Concatenate all individual DataFrames
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'{path}/MLE_summary_{num_simulations}simul_{n_bins}bins_{max_mut}MaxMut_{distrib}.csv', index=False)

########################################################################################################################
