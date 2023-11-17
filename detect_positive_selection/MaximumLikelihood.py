#!/usr/bin/env python
# coding=utf-8

import numpy as np
from scipy import stats
from scipy.optimize import minimize
import matplotlib as mat
import matplotlib.pyplot as plt
from scipy.stats import chi2
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D


def proba_fixation(delta, params):
    if len(params) == 0:
        return 1
    else:
        std_dev = params[0]
        mean = params[1] if len(params) == 2 else 0
        w_mutant = stats.norm.pdf(delta, mean, std_dev)
        w_ancestral = stats.norm.pdf(0, mean, std_dev)
        s = (w_mutant - w_ancestral)/w_ancestral
        print(s)
        if not np.isfinite(s):
            return 0
        pfix = s/(1-np.exp(-s)) if abs(s) > 1.e-4 else 1+(s/2)
        return pfix


def proba_substitution(params, hist_mutations):
    n_bins = len(hist_mutations[0])
    output_array = np.zeros(n_bins)
    for b in range(n_bins):
        p = hist_mutations[0][b]
        min_bin = hist_mutations[1][b]
        max_bin = hist_mutations[1][b+1]
        delta = (max_bin - min_bin) / 2
        output_array[b] = p * proba_fixation(delta, params)
    output_array /= np.sum(output_array)
    return output_array, hist_mutations[1]


def find_bins(hist, value):
    n_bins = len(hist[0])
    for b in range(n_bins):
        min_bin = hist[1][b]
        max_bin = hist[1][b + 1]
        if min_bin < value < max_bin:
            return hist[0][b]
    return 0.0


def loglikelihood(deltas, params, hist_mutations):
    hist_subs = proba_substitution(params, hist_mutations)
    models_lk = []
    for delta in deltas:
        models_lk.append(find_bins(hist_subs, delta))
    if np.sum(models_lk) == 0:
        return np.infty
    ln_lk = -np.sum(np.log(models_lk))
    return ln_lk


# Plot the results for each sequence
def general_plot(all_deltas, obs, svm_distribution, purif, pos, scenario, test, sub_ax1, sub_ax2):
    sub_ax1.hist(all_deltas, bins=50, density=True, alpha=0.7, label="All mutations")
    sub_ax1.plot(sorted(all_deltas), svm_distribution(sorted(all_deltas)))
    for obs_value in obs:
        sub_ax1.axvline(obs_value, color='red', linestyle='--')
    sub_ax1.set_xlabel("deltaSVM")
    sub_ax1.set_ylabel("Density")

    legend_elements = [mat.lines.Line2D([0], [0], color='blue', linestyle='-', linewidth=0.1,
                                        label=f'All: m={round(np.mean(all_deltas),2)}, std={round(np.std(all_deltas),2)}'),
                       mat.lines.Line2D([0], [0], color='red', linestyle='--', linewidth=0.1,
                                        label=f'Sub: m={round(np.mean(obs),2)}, std={round(np.std(obs),2)}')]
    sub_ax1.legend(handles=legend_elements)
    sub_ax1.set_title(f'Simulation for a {scenario} case')

    # Plot of observed Sub and fitted distribution
    x_values = np.linspace(min(all_deltas), max(all_deltas), 1000)
    density_uniform_constant = 1 / (max(all_deltas) - min(all_deltas))
    density_uniform = np.full_like(x_values, density_uniform_constant)
    density_purif = stats.norm.pdf(x_values, 0, purif.x[0])
    density_pos = stats.norm.pdf(x_values, pos.x[1], pos.x[0])

    sub_ax2.plot(x_values, density_uniform, color="black", label="Random Drift")
    sub_ax2.plot(x_values, density_purif, color="red", label=f'Purif: m=0, std={round(purif.x[0], 2)}')
    sub_ax2.plot(x_values, density_pos, color="green", label=f'Pos: m={round(pos.x[1], 2)}, std={round(pos.x[0], 2)}')
    for obs_value in obs:
        sub_ax2.axvline(obs_value, color='red', linestyle='-', linewidth=0.01)

    sub_ax2.set_xlabel("deltaSVM")
    sub_ax2.set_ylabel("Density")
    sub_ax2.legend()
    sub_ax2.set_title(f'Maximum Likelihood Estimation: {test}')


# Plot the minimization process for each sequence
def plot_model(obs, model_params, ax, model_type="Stabilizing"):
    theta_range = np.linspace(0.000001, 10, 200)  # Std range

    if model_type == "Stabilizing":
        ll_values = [loglikelihood(obs, [theta], hist_SVM) for theta in theta_range]
        ax.plot(theta_range, ll_values, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], loglikelihood(obs, model_params, hist_SVM), color='red', marker='o', label="Minimized Point")
        ax.set_xlabel("Parameter: Standard Deviation")
        ax.set_ylabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        ax.legend()

    elif model_type == "Positive":
        theta_range = np.linspace(0.000001, 2, 70)  # Std range
        theta2_range = np.linspace(-5, 7, 70)  # Mean range

        theta1_values, theta2_values = np.meshgrid(theta_range, theta2_range)
        ll_values = np.array([[loglikelihood(obs, [theta1, theta2], hist_SVM) for theta1 in theta_range] for theta2 in theta2_range])

        # Plotting the log-likelihood surface
        ax.plot_surface(theta1_values, theta2_values, ll_values, cmap='viridis', alpha=0.8, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], model_params[1], loglikelihood(obs, model_params, hist_SVM), color='red', marker='o', label="Minimized Point")
        ax.invert_xaxis()
        ax.set_xlabel("Parameter: Variance")
        ax.set_ylabel("Parameter: Mean")
        ax.set_zlabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        #ax.legend()

########################################################################################################################
# Datas
path = "/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/"
num_simulations = 2
dfs = []
for Nsim in range(num_simulations):
    # SVM score distribution: affinity of all possible deltas for a sequence
    All_SVM = np.random.normal(-1, 3, 2500)
    gaussian_mutation = stats.gaussian_kde(All_SVM)    # fit a kernel distribution on our data
    hist_SVM = np.histogram(All_SVM, bins=100)

    # Observed deltaSVM: change in affinity of substitutions for a sequence
    Nsub = np.random.randint(5, 25)
    Simul_Obs_SVM = {"Random": np.random.normal(-1, 3, Nsub),
                     "Stabilizing": np.random.normal(0, 1, Nsub),
                     "Positive": np.random.normal(5, 1, Nsub)}

    for simul in Simul_Obs_SVM.keys():
        print(f'Simulation_{Nsim}: scenario {simul}')
        obs_svm = Simul_Obs_SVM[simul]

        # Null model: no param.
        LL_neutral = loglikelihood(obs_svm, [], hist_SVM)

        # Stabilizing selection: fixed mean to 0, minimize variance
        model_purif = minimize(lambda theta: loglikelihood(obs_svm, theta, hist_SVM), np.array([1]),
                               bounds=[(0, np.infty)], method="SLSQP")
        LL_purif = model_purif.fun
        print(f"Stabilizing model:{model_purif}")

        # Positive selection: minimize mean AND variance
        model_pos = minimize(lambda theta: loglikelihood(obs_svm, theta, hist_SVM), np.array([1, 0]),
                             bounds=[(0, np.infty), (-np.infty, np.infty)], method="SLSQP")
        LL_pos = model_pos.fun
        print(f"Positive model:{model_pos}")

        # Likelihood ratio test:
        LRT_null_purif = 2 * (LL_neutral-LL_purif)
        LRT_null_pos = 2 * (LL_neutral - LL_pos)
        LRT_purif_pos = 2 * (LL_purif - LL_pos)

        p_value_null_purif = chi2.sf(LRT_null_purif, 1)
        p_value_null_pos = chi2.sf(LRT_null_pos, 2)
        p_value_purif_pos = chi2.sf(LRT_purif_pos, 1)

        if p_value_null_purif < 0.05:
            conclusion = "Stabilizing model"
            if p_value_purif_pos < 0.05:
                conclusion = "Positive model"
        else:
            if p_value_null_pos < 0.05:
                conclusion = "Positive model"
            else:
                conclusion = "Null model"

        print("Conclusion:", conclusion)

        # Create a DataFrame for each simulation
        df = pd.DataFrame({"Nmut": [len(obs_svm)], "MeanObs": [np.mean(obs_svm)], "VarObs": [np.var(obs_svm)],
                           "Scenario": [simul], "VarPurif:": [model_purif.x[0]], "VarPos:": [model_pos.x[0]],
                           "MeanPos:": [model_pos.x[1]], "NiterPurif:": [model_purif.nit], "NiterPos:": [model_pos.nit],
                           "LL_neutral": [LL_neutral], "LL_purif": [LL_purif], "LL_pos": [LL_pos],
                           "Conclusion": [conclusion]})

        dfs.append(df)


        # with PdfPages(f'{path}/MLE_summary_{simul}_{Nsim}_v2.pdf') as pdf:
        #     fig, axes = plt.subplots(2, 2, figsize=(14, 14))
        #     axes[1, 1].axis('off')
        #     axes[1, 1] = fig.add_subplot(224, projection='3d')
        #     general_plot(All_SVM, obs_svm, gaussian_mutation, model_purif, model_pos, simul, conclusion, axes[0,0], axes[0, 1])
        #     plot_model(obs_svm, model_purif.x, axes[1, 0], model_type="Stabilizing")
        #     plot_model(obs_svm, model_pos.x, axes[1, 1], model_type="Positive")
        #     pdf.savefig()
        #     plt.close(fig)


# Concatenate all individual DataFrames
result_df = pd.concat(dfs, ignore_index=True)

# Write the DataFrame to a CSV file
#result_df.to_csv('/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/simulation_results.csv', index=False)

########################################################################################################################





