#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2


# Get the quantiles of the SVM distribution of each side
def get_svm_quantiles(all_svm, obs_svm, n_quant=25):
    neg_svm = [x for x in all_svm if x <= 0]
    pos_svm = [x for x in all_svm if x > 0]

    obs_bins = []
    # Ensure that all substitutions don't fall into the same quantile
    while True:
        # split the distribution into Nb quantiles on each side of the distribution
        neg_quant, neg_bins = pd.qcut(neg_svm, q=n_quant, retbins=True)
        pos_quant, pos_bins = pd.qcut(pos_svm, q=n_quant, retbins=True)
        bins_values = list(neg_bins) + list(pos_bins)[1:]  # merge bins and remove the first of pos_bins

        # Get the quantile of the observed values
        obs_bins = np.searchsorted(bins_values, obs_svm, side='left') - 1
        if len(set(obs_bins)) >= 2:
            break
        n_quant += 5

    n_bins = n_quant * 2
    # Get the probability of each quantile
    quant_count = neg_quant.value_counts().to_list() + pos_quant.value_counts().to_list()
    quant_proba = quant_count / np.sum(quant_count)
    obs_quant = [1 / (2 * n_bins) + b / n_bins for b in obs_bins]
    return quant_proba, obs_bins, obs_quant


# Selection coefficient from delta's quantile and model's parameters
def coeff_selection(quant_delta, params):
    assert 0 < quant_delta < 1
    alpha = params[0] if len(params) > 0 else 1.0
    beta = params[1] if len(params) == 2 else alpha
    w_mutant = stats.beta.pdf(quant_delta, a=alpha, b=beta)
    w_ancestral = stats.beta.pdf(0.5, a=alpha, b=beta)
    if w_ancestral < 1.e-10:
        print(f"q:{quant_delta:.2f}, w_m:{w_mutant:.2g}, w_a:{w_ancestral:.2g}, alpha:{alpha:.2f}, beta:{beta:.2f}")
    if w_mutant < 1.e-10 or w_ancestral < 1.e-10:
        return -np.infty
    s = np.log(w_mutant / w_ancestral)
    return s


# Scaled fixation probability from selection coefficient
def proba_fixation(s):
    if s == -np.infty:
        return 0.0
    elif abs(s) < 1.e-4:
        return 1.0 + s / 2
    elif s > 10:
        return s
    elif s < -10:
        return - s * np.exp(s)
    else:
        return s / (1 - np.exp(-s))


# Get probability of substitution for each quantile: P(Mut) * P(Fix)
def proba_substitution(params, mutations_proba):
    n_bins = len(mutations_proba)
    output_array = np.zeros(n_bins)
    delta = 1 / (2 * n_bins)
    for quant in range(n_bins):
        # Probability of mutation
        proba_mut = mutations_proba[quant]

        # Probability of fixation
        quant_val = quant / n_bins + delta  # quantile needs to be between 0<quant<1
        s = coeff_selection(quant_val, params)
        proba_fix = proba_fixation(s)

        output_array[quant] = proba_mut * proba_fix

    # Scaled output
    sum_output = np.sum(output_array)
    if sum_output == 0:
        return output_array
    else:
        return output_array / sum_output


def loglikelihood(obs_svm, params, all_svm, lambda_rate=0.0):
    if len(params) == 1 and params[0] < 1e-10:
        return -np.infty
    if len(params) == 2 and (params[0] < 1e-10 or params[1] < 1e-10):
        return -np.infty
    mutations_proba, obs_bins, obs_quant = get_svm_quantiles(all_svm, obs_svm)
    subs_proba = proba_substitution(params, mutations_proba)
    models_lk = [subs_proba[i] if i >= 0 else subs_proba[0] for i in obs_bins]

    if min(models_lk) <= 0.0:  # Avoid log(0)
        return -np.infty
    lnl = np.sum(np.log(models_lk))

    # Add a penalty params far from 1.0
    if lambda_rate > 0.0 and len(params) == 1:
        lnl -= lambda_rate * abs(np.log(params[0]))
    if lambda_rate > 0.0 and len(params) == 2:
        lnl -= lambda_rate * (abs(np.log(params[0])) + abs(np.log(params[1])) + abs(np.log(params[0] / params[1]))) / 3

    return lnl


def AIC_weight(ll_neutral, ll_purif, ll_pos, df_purif, df_pos):
    aic_neutral = -2 * ll_neutral
    aic_purif = -2 * ll_purif + 2 * df_purif
    aic_pos = -2 * ll_pos + 2 * df_pos
    max_aic = min(aic_neutral, aic_purif, aic_pos)
    weights = np.exp(-0.5 * (np.array([aic_neutral, aic_purif, aic_pos]) - max_aic))
    return weights / np.sum(weights)


def conclusion_purif(alpha):
    return "Stabilizing (m=0.5)" if alpha > 1.0 else "Disruptive (m=0.5)"


def conclusion_pos(alpha, beta):
    if alpha > 1.0 and beta > 1.0:
        return f"Stabilizing (m={alpha/(alpha+beta):.2f})"
    if alpha > 1.0 > beta:
        return "Directional (+)"
    if alpha < 1.0 < beta:
        return "Directional (-)"
    else:
        return f"Disruptive (m={alpha/(alpha+beta):.2f})"


def run_estimations(all_svm, obs_svm, alpha_threshold=0.05, lambda_rate=0.05, verbose=False, ID=""):
    # Null model: no param.
    ll_neutral = loglikelihood(obs_svm, [], all_svm)

    # Stabilizing selection: alpha = beta
    bounds = [(0.0, np.inf)]
    model_purif = minimize(lambda theta: -loglikelihood(obs_svm, theta, all_svm, lambda_rate), np.array([1.0]),
                           bounds=bounds, method="Nelder-Mead")
    ll_purif = loglikelihood(obs_svm, model_purif.x, all_svm)

    # Positive selection
    bounds = [(0.0, np.inf), (0.0, np.inf)]
    initial_guess = np.array([1.0, 1.0])
    model_pos = minimize(lambda theta: -loglikelihood(obs_svm, theta, all_svm, lambda_rate), initial_guess,
                         bounds=bounds, method="Nelder-Mead")
    ll_pos = loglikelihood(obs_svm, model_pos.x, all_svm)

    # Likelihood ratio test:
    lrt_null_purif = -2 * (ll_neutral - ll_purif)
    lrt_purif_pos = -2 * (ll_purif - ll_pos)

    p_value_null_purif = chi2.sf(lrt_null_purif, 1)
    p_value_purif_pos = chi2.sf(lrt_purif_pos, 1)
    p_value_null_pos = chi2.sf(lrt_purif_pos, 2)
    conclusion = "Neutral"
    if p_value_null_purif < alpha_threshold:
        conclusion = conclusion_purif(model_purif.x[0])
    if p_value_purif_pos < p_value_null_purif and p_value_purif_pos < alpha_threshold:
        conclusion = conclusion_pos(model_pos.x[0], model_pos.x[1])

    df_purif = 1 if (lambda_rate == 0.0) else int(abs(model_purif.x[0] - 1.0) > 1.e-10)
    df_pos = 2 if (lambda_rate == 0.0) else (int(abs(model_pos.x[0] - 1.0) > 1.e-10) +
                                             int(abs(model_pos.x[1] - 1.0) > 1.e-10))
    weight_neutral, weight_purif, weight_pos = AIC_weight(ll_neutral, ll_purif, ll_pos, df_purif, df_pos)
    conclusion_weighted = "None"
    if weight_neutral > max(weight_purif, weight_pos):
        conclusion_weighted = "Neutral"
    if weight_purif > max(weight_neutral, weight_pos):
        conclusion_weighted = conclusion_purif(model_purif.x[0])
    if weight_pos > max(weight_neutral, weight_purif):
        conclusion_weighted = conclusion_pos(model_pos.x[0], model_pos.x[1])

    if verbose:
        print(f"  LnL neutral    : {ll_neutral:.2f}")
        print(f"  LnL stabilizing: {ll_purif:.4f}, alpha: {model_purif.x[0]:.3g}, beta: {model_purif.x[0]:.3g}")
        print(f"  LnL positive   : {ll_pos:.4f}, alpha: {model_pos.x[0]:.3g}, beta: {model_pos.x[1]:.3g}")
        print(f"  LRT null vs purif: {lrt_null_purif:.2f}, p-value: {p_value_null_purif:.2g}")
        print(f"  LRT purif vs pos : {lrt_purif_pos:.2f}, p-value: {p_value_purif_pos:.2g}")
        print("  Conclusion:", conclusion)

    # Create a DataFrame
    result = pd.DataFrame(
        {"ID": [ID], "Nmut": [len(obs_svm)], "SumObs": [np.sum(obs_svm)], "MeanObs": [np.mean(obs_svm)],
         "VarObs": [np.var(obs_svm)], "MedSVM": [np.median(all_svm)],
         "MinSVM": [np.min(all_svm)], "MaxSVM": [np.max(all_svm)],
         "AlphaPurif": [model_purif.x[0]], "AlphaPos": [model_pos.x[0]], "BetaPos": [model_pos.x[1]],
         "NdfPurif": [df_purif], "NdfPos": [df_pos],
         "NiterPurif": [model_purif.nit], "NiterPos": [model_pos.nit],
         "LL_neutral": [ll_neutral], "LL_purif": [ll_purif], "LL_pos": [ll_pos],
         "LRT_null_purif": [lrt_null_purif], "LRT_purif_pos": [lrt_purif_pos],
         "p_value_null_purif": [p_value_null_purif], "p_value_purif_pos": [p_value_purif_pos],
         "p_value_null_pos": [p_value_null_pos], "Conclusion": [conclusion],
         "w_neutral": [weight_neutral], "w_purif": [weight_purif], "w_pos": [weight_pos],
         "Conclusion_weighted": [conclusion_weighted]})

    return result, [model_purif, model_pos, bounds]


# Plot the results for each sequence
def general_plot(all_deltas, obs, svm_distribution, purif, pos, scenario, test, sub_ax1, sub_ax2):
    sub_ax1.hist(all_deltas, bins=50, density=True, alpha=0.7, label="All mutations")
    sub_ax1.plot(sorted(all_deltas), svm_distribution(sorted(all_deltas)))
    for obs_value in obs:
        sub_ax1.axvline(obs_value, color='red', linestyle='--')
    sub_ax1.axvline(0.0, color='black', linestyle='-')
    sub_ax1.set_xlabel("deltaSVM")
    sub_ax1.set_ylabel("Density")

    legend_elements = [Line2D([0], [0], color='blue', linestyle='-', linewidth=0.1,
                              label=f'All: m={round(np.mean(all_deltas), 2)}, std={round(np.std(all_deltas), 2)}'),
                       Line2D([0], [0], color='red', linestyle='--', linewidth=0.1,
                              label=f'Sub: m={round(np.mean(obs), 2)}, std={round(np.std(obs), 2)}')]
    sub_ax1.legend(handles=legend_elements)
    sub_ax1.set_title(f'{scenario}')

    # Plot of observed Sub and fitted distribution
    mutations_proba, obs_bins, obs_quant = get_svm_quantiles(all_deltas, obs)
    neu_subs_proba = proba_substitution([], mutations_proba)
    purif_subs_proba = proba_substitution(purif.x, mutations_proba)
    pos_subs_proba = proba_substitution(pos.x, mutations_proba)

    # Bar plot neutral subs
    sub_ax2.plot(np.arange(len(neu_subs_proba)), neu_subs_proba, color='blue', alpha=0.5,
                 label=f'Neutral: $\\alpha=\\beta=1.0$')
    sub_ax2.plot(np.arange(len(neu_subs_proba)), purif_subs_proba, color="red",
                 label=f'Purif: $\\alpha=\\beta={purif.x[0]:.2g}$')
    sub_ax2.plot(np.arange(len(neu_subs_proba)), pos_subs_proba, color="green",
                 label=f'Pos: $\\alpha={pos.x[0]:.2g}, \\beta={pos.x[1]:.2g}$')

    for obs_value in obs_bins:
        sub_ax2.axvline(obs_value, color='red', linestyle='-', linewidth=0.01)
    sub_ax2.axvline(len(neu_subs_proba) // 2, color='black', linestyle='-', linewidth=0.05)

    sub_ax2.set_xlabel("deltaSVM")
    sub_ax2.set_ylabel("Density")
    sub_ax2.legend()
    sub_ax2.set_title(f'Maximum Likelihood Estimation: {test}')


# Plot the minimization process for each sequence
def plot_model(obs, all_svm, model_params, ax, model_type="Stabilizing", bounds=None, lambda_rate=0.0):
    if model_type == "Stabilizing":
        alpha_min = max([bounds[0][0], model_params[0] * 0.5])
        alpha_max = min([bounds[0][1], model_params[0] * 1.5])
        alpha_range = np.linspace(alpha_min, alpha_max, 50)  # Std range
        ll_values = [loglikelihood(obs, [a], all_svm, lambda_rate) for a in alpha_range]
        ax.plot(alpha_range, ll_values, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], loglikelihood(obs, model_params, all_svm, lambda_rate), color='red', marker='o',
                   label="Maximum likelihood")
        ax.set_xlabel("Parameter: $\\alpha=\\beta$")
        ax.set_ylabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        ax.legend()

    elif model_type == "Directional":
        alpha_min = max([bounds[0][0], model_params[0] * 0.5])
        alpha_max = min([bounds[0][1], model_params[0] * 1.5])
        beta_min = max([bounds[1][0], model_params[1] * 0.5])
        beta_max = min([bounds[1][1], model_params[1] * 1.5])
        alpha_range = np.linspace(alpha_min, alpha_max, 10)
        beta_range = np.linspace(beta_min, beta_max, 10)

        alpha_values, beta_values = np.meshgrid(alpha_range, beta_range)
        ll_values = np.array(
            [[loglikelihood(obs, [a, b], all_svm, lambda_rate) for a in alpha_range] for b in beta_range])

        # Plotting the log-likelihood surface
        ax.plot_surface(alpha_values, beta_values, ll_values, cmap='viridis', alpha=0.8,
                        label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], model_params[1], loglikelihood(obs, model_params, all_svm, lambda_rate),
                   color='red',
                   marker='o', label="Maximum likelihood")
        ax.invert_xaxis()
        ax.set_xlabel("Parameter: $\\alpha$")
        ax.set_ylabel("Parameter: $\\beta$")
        ax.set_zlabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")


################################# !!!!!!! Temporary to test manually !!!!!!! ###########################################
DeltaSVM = pd.read_csv("ancestral_all_possible_deltaSVM.txt", sep='\t', header=0)
AllObsSVM = pd.read_csv("ancestral_to_observed_deltaSVM.txt", sep='\t', header=None, names=range(150 + 4))
AllObsSVM.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(AllObsSVM.columns[4:])

# ID = "chr13:86947569:86947706_13:86947569:86947706:Interval_6751"
for penalized_rate in [0.0, 0.05, 0.1, 0.2]:
    os.makedirs(f"plots_{penalized_rate}", exist_ok=True)
    result_list = []
    for ID in AllObsSVM['ID']:
        all_svm_row = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0]
        obs_svm_row = AllObsSVM.loc[AllObsSVM['ID'] == ID, 4:].iloc[0]
        all_svm = all_svm_row.dropna().values.tolist()
        obs_svm = obs_svm_row.dropna().values.tolist()
        print(f"Obs SVM:{obs_svm}")
        results, models = run_estimations(all_svm, obs_svm, 0.05, penalized_rate, ID=ID)
        print(results)
        fig, axes = plt.subplots(2, 2, figsize=(14, 14))
        axes[1, 1].axis('off')
        axes[1, 1] = fig.add_subplot(224, projection='3d')
        gaussian_mutation = stats.gaussian_kde(all_svm)
        general_plot(all_svm, obs_svm, gaussian_mutation, models[0], models[1], ID, results["Conclusion_weighted"][0],
                     axes[0, 0], axes[0, 1])

        plot_model(obs_svm, all_svm, models[0].x, axes[1, 0], model_type="Stabilizing", bounds=models[2],
                   lambda_rate=penalized_rate)
        plot_model(obs_svm, all_svm, models[1].x, axes[1, 1], model_type="Directional", bounds=models[2],
                   lambda_rate=penalized_rate)
        plt.savefig(f'plots_{penalized_rate}/MLE_summary_{ID}.pdf')
        plt.close(fig)
        plt.clf()
        result_list.append(results)

    final_results = pd.concat(result_list)
    final_results.to_csv(f"MLE_results_{penalized_rate}.csv", index=False)
    ########################################################################################################################
