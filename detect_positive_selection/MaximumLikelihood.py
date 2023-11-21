#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
from scipy.optimize import minimize
import matplotlib as mat
import matplotlib.pyplot as plt
from scipy.stats import chi2
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

# Seed numpy random generator
np.random.seed(1234)


def proba_fixation(delta, params):
    if len(params) == 0:
        return 1
    else:
        std_dev = params[0]
        if std_dev == 0:
            return 0.0
        mean = params[1] if len(params) == 2 else 0
        w_mutant = stats.norm.pdf(delta, loc=mean, scale=std_dev)
        w_ancestral = stats.norm.pdf(mean, loc=mean, scale=std_dev)
        if w_mutant == 0.0:
            return 0.0
        s = np.log(w_mutant / w_ancestral)
        if abs(s) < 1.e-4:
            return 1.0 + s / 2
        elif s > 10:
            return s
        elif s < -10:
            return - s * np.exp(s)
        else:
            return s / (1 - np.exp(-s))


def proba_substitution(params, mutations_proba, bins_values):
    n_bins = len(mutations_proba)
    output_array = np.zeros(n_bins)
    for b in range(n_bins):
        p = mutations_proba[b]
        min_bin = bins_values[b]
        max_bin = bins_values[b + 1]
        mean_bin = (max_bin + min_bin) / 2
        output_array[b] = p * proba_fixation(mean_bin, params)
    sum_output = np.sum(output_array)
    if sum_output == 0:
        return output_array
    else:
        return output_array / sum_output


def loglikelihood(deltas, params, hist_mutations):
    bins_values = hist_mutations[1]
    mutations_proba = hist_mutations[0] / np.sum(hist_mutations[0])
    subs_proba = proba_substitution(params, mutations_proba, bins_values)

    digitized_deltas = np.searchsorted(bins_values, deltas, side='left') - 1
    models_lk = [subs_proba[i] if i >= 0 else subs_proba[0] for i in digitized_deltas]

    if min(models_lk) <= 0.0:
        return -np.infty
    return np.sum(np.log(models_lk))


# Plot the results for each sequence
def general_plot(all_deltas, obs, svm_distribution, purif, pos, scenario, test, sub_ax1, sub_ax2):
    sub_ax1.hist(all_deltas, bins=50, density=True, alpha=0.7, label="All mutations")
    sub_ax1.plot(sorted(all_deltas), svm_distribution(sorted(all_deltas)))
    for obs_value in obs:
        sub_ax1.axvline(obs_value, color='red', linestyle='--')
    sub_ax1.set_xlabel("deltaSVM")
    sub_ax1.set_ylabel("Density")

    legend_elements = [mat.lines.Line2D([0], [0], color='blue', linestyle='-', linewidth=0.1,
                                        label=f'All: m={round(np.mean(all_deltas), 2)}, std={round(np.std(all_deltas), 2)}'),
                       mat.lines.Line2D([0], [0], color='red', linestyle='--', linewidth=0.1,
                                        label=f'Sub: m={round(np.mean(obs), 2)}, std={round(np.std(obs), 2)}')]
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
def plot_model(obs, model_params, ax, model_type="Stabilizing", bounds=None):
    if model_type == "Stabilizing":
        k = 0.25
        std_min = max([bounds[0][0], model_params[0] * k])
        std_max = min([bounds[0][1], model_params[0] * (2 - k)])
        std_range = np.linspace(std_min, std_max, 100)  # Std range
        ll_values = [loglikelihood(obs, [std], hist_SVM) for std in std_range]
        ax.plot(std_range, ll_values, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], loglikelihood(obs, model_params, hist_SVM), color='red', marker='o',
                   label="Minimized Point")
        ax.set_xlabel("Parameter: Standard Deviation")
        ax.set_ylabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        ax.legend()

    elif model_type == "Positive":
        k = 0.25
        std_min = max([bounds[0][0], model_params[0] * k])
        std_max = min([bounds[0][1], model_params[0] * (2 - k)])
        mean_min = max([bounds[1][0], model_params[1] - 3])
        mean_max = min([bounds[1][1], model_params[1] + 3])
        std_range = np.linspace(std_min, std_max, 30)
        mean_range = np.linspace(mean_min, mean_max, 30)

        std_values, mean_values = np.meshgrid(std_range, mean_range)
        ll_values = np.array(
            [[loglikelihood(obs, [std, mean], hist_SVM) for std in std_range] for mean in mean_range])

        # Plotting the log-likelihood surface
        ax.plot_surface(std_values, mean_values, ll_values, cmap='viridis', alpha=0.8,
                        label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], model_params[1], loglikelihood(obs, model_params, hist_SVM), color='red',
                   marker='o', label="Minimized Point")
        ax.invert_xaxis()
        ax.set_xlabel("Parameter: Variance")
        ax.set_ylabel("Parameter: Mean")
        ax.set_zlabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        # ax.legend()


########################################################################################################################
# Datas
path = "results/"
os.makedirs(path, exist_ok=True)
num_simulations = 1000
n_bins = 30
no_plots = True
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
    Nsub = np.random.randint(5, 20)
    print(f'Number of substitutions: {Nsub}')
    Simul_Obs_SVM = {"Random": np.random.choice(All_SVM, Nsub, replace=False),
                     "Stabilizing": np.random.choice(All_SVM, Nsub, replace=False, p=stab_weights),
                     "Positive": np.random.choice(All_SVM, Nsub, replace=False, p=pos_weights)}

    for simul in Simul_Obs_SVM.keys():
        print(f'Simulation {simul} {Nsim}')
        obs_svm = Simul_Obs_SVM[simul]

        # Null model: no param.
        LL_neutral = loglikelihood(obs_svm, [], hist_SVM)
        print(f"  LnL neutral    : {LL_neutral:.2f}")

        std_max = np.std((min(All_SVM), max(All_SVM)))
        std_min = 1e-1
        # Stabilizing selection: fixed mean to 0, minimize variance
        model_purif = minimize(lambda theta: -loglikelihood(obs_svm, theta, hist_SVM), np.array([std_max / 2]),
                               bounds=[(std_min, std_max)], method="L-BFGS-B")
        LL_purif = -model_purif.fun
        print(f"  LnL stabilizing: {LL_purif:.2f}, std: {model_purif.x[0]:.2f}")

        optimal_min = hist_SVM[1][1]
        optimal_max = hist_SVM[1][-2]
        bounds = [(std_min, std_max), (optimal_min, optimal_max)]
        # Positive selection: minimize mean AND variance
        model_pos = minimize(lambda theta: -loglikelihood(obs_svm, theta, hist_SVM), np.array([model_purif.x[0], 0]),
                             bounds=bounds, method="Powell")
        LL_pos = -model_pos.fun
        print(f"  LnL positive   : {LL_pos:.2f}, std: {model_pos.x[0]:.2f}, mean: {model_pos.x[1]:.2f}")

        # Likelihood ratio test:
        LRT_null_purif = -2 * (LL_neutral - LL_purif)
        LRT_null_pos = -2 * (LL_neutral - LL_pos)
        LRT_purif_pos = -2 * (LL_purif - LL_pos)

        p_value_null_purif = chi2.sf(LRT_null_purif, 1)
        p_value_null_pos = chi2.sf(LRT_null_pos, 2)
        p_value_purif_pos = chi2.sf(LRT_purif_pos, 1)
        print(f"  LRT null vs purif: {LRT_null_purif:.2f}, p-value: {p_value_null_purif:.2g}")
        print(f"  LRT null vs pos  : {LRT_null_pos:.2f}, p-value: {p_value_null_pos:.2g}")
        print(f"  LRT purif vs pos : {LRT_purif_pos:.2f}, p-value: {p_value_purif_pos:.2g}")

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

        print("  Conclusion:", conclusion)

        # Create a DataFrame for each simulation
        df = pd.DataFrame({"Nmut": [len(obs_svm)], "MeanObs": [np.mean(obs_svm)], "VarObs": [np.var(obs_svm)],
                           "Scenario": [simul], "VarPurif:": [model_purif.x[0]], "VarPos:": [model_pos.x[0]],
                           "MeanPos:": [model_pos.x[1]], "NiterPurif:": [model_purif.nit], "NiterPos:": [model_pos.nit],
                           "LL_neutral": [LL_neutral], "LL_purif": [LL_purif], "LL_pos": [LL_pos],
                           "Conclusion": [conclusion]})

        dfs.append(df)
        if no_plots:
            continue
        with PdfPages(f'{path}/MLE_summary_{simul}_{Nsim}.pdf') as pdf:
            fig, axes = plt.subplots(2, 2, figsize=(14, 14))
            axes[1, 1].axis('off')
            axes[1, 1] = fig.add_subplot(224, projection='3d')
            general_plot(All_SVM, obs_svm, gaussian_mutation, model_purif, model_pos, simul, conclusion, axes[0, 0],
                         axes[0, 1])
            plot_model(obs_svm, model_purif.x, axes[1, 0], model_type="Stabilizing", bounds=bounds)
            plot_model(obs_svm, model_pos.x, axes[1, 1], model_type="Positive", bounds=bounds)
            pdf.savefig()
            plt.close(fig)
            plt.clf()

# Concatenate all individual DataFrames
result_df = pd.concat(dfs, ignore_index=True)

# Write the DataFrame to a CSV file
result_df.to_csv(f'{path}/MLE_summary.csv', index=False)

########################################################################################################################