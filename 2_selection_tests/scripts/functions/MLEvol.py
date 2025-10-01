#!/usr/bin/env python
# coding=utf-8
import numpy as np
from scipy import stats
import matplotlib as mat
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2
import numba as nb
np.random.seed(1234)


def sum_nuc_probs(nuc_changes, sub_mat_proba):
    # Get the nucleotide probabilities for each substitution
    sum_rate = 0
    for i in range(len(nuc_changes)):
        pos, source, target = nuc_changes[i].split(":")
        sum_rate += sub_mat_proba[source][target]
    return sum_rate


# Get SVM ids (pos:nuc:sub)
def get_svm_ids(sequence):
    svm_ids = []
    for pos_id in sequence:
        pos = pos_id.split(":")[0]
        nuc = pos_id.split(":")[1]
        for sub in ["A", "T", "C", "G"]:
            if nuc != sub:
                svm_ids.append(pos + ":" + nuc + ":" + sub)
    return svm_ids

# Get the quantiles of the SVM distribution of each side
def get_svm_quantiles(all_svm, obs_svm, all_svm_ids, sub_mat_proba, n_quant=50):
    neg_svm = [x for x in all_svm if x <= 0]
    pos_svm = [x for x in all_svm if x > 0]

    obs_bins = []
    # Ensure that all substitutions don't fall into the same quantile
    while len(set(obs_bins)) < 2:
        # split the distribution into Nb quantiles/2 on each side of the distribution
        neg_quant, neg_bins = pd.qcut(neg_svm, q=int(n_quant/2), retbins=True, duplicates='drop')
        pos_quant, pos_bins = pd.qcut(pos_svm, q=int(n_quant/2), retbins=True, duplicates='drop')
        bins_values = list(neg_bins) + list(pos_bins)[1:]  # merge bins and remove the first of pos_bins

        # Get the quantiles of the observed values
        obs_bins = np.searchsorted(bins_values, obs_svm, side='left') - 1
        n_quant += 5
        if n_quant > 200:
            break

    # Get the probability of each quantile
    quant_count = neg_quant.value_counts().to_list() + pos_quant.value_counts().to_list()
    # Sort the SVM values and keep the nuc_changes in the same order
    sorted_svm, sorted_changes = zip(*sorted(zip(all_svm, all_svm_ids)))
    min_quant = 0
    for i in range(len(quant_count)):
        # extract the SVM changes for each quantile
        max_quant = min_quant + quant_count[i]
        quant_changes = sorted_changes[min_quant:max_quant]
        assert len(quant_changes) == quant_count[i]
        # Get the nucleotide probabilities for each substitution
        quant_count[i] = sum_nuc_probs(quant_changes, sub_mat_proba)
        min_quant = max_quant

    quant_proba = quant_count / np.sum(quant_count)

    # Scaling the bins values to be between 0 and 1
    delta = 1 / (2 * n_quant)
    scaled_bins = [(i / n_quant) + delta for i in range(len(quant_proba))]

    return quant_proba, obs_bins, scaled_bins


def get_svm_hist(all_svm, obs_svm, n_bin=50):
    obs_bins = []
    # Ensure that all substitutions don't fall into the same bin
    while len(set(obs_bins)) < 2:
        hist_svm = np.histogram(all_svm, bins=n_bin)
        bins_values = hist_svm[1]

        # Get the bin of the observed values
        obs_bins = np.searchsorted(bins_values, obs_svm, side='left') - 1
        n_bin += 5
        if n_bin > 200:
            break

    # Get the probability of each bin
    bin_proba = hist_svm[0] / np.sum(hist_svm[0])

    # Scaling the bins values to be centered on 0.5
    scaled_bins = np.zeros(len(bin_proba))
    delta_bounds = [min(bins_values), max(bins_values)]
    for b in range(len(bin_proba)):
        min_bin = bins_values[b]
        max_bin = bins_values[b + 1]
        mean_bin = (max_bin + min_bin) / 2
        max_delta = max(np.abs(delta_bounds))
        scaled_bins[b] = (mean_bin + max_delta) / (2 * max_delta)

    return bin_proba, obs_bins, scaled_bins


def get_obs_index(obs_svm, sorted_svm, flag=""):
    obs_svm = np.atleast_1d(obs_svm)
    obs_index = np.atleast_1d(np.searchsorted(sorted_svm, obs_svm, side='left'))
    for i in range(len(obs_index)):
        if obs_index[i] < 0:
            obs_index[i] = 0
        if flag != "alt_ancestral":
            assert sorted_svm[obs_index[i]] == obs_svm[i]
    assert len(obs_index) == len(obs_svm)
    assert min(obs_index) >= 0
    assert max(obs_index) < len(sorted_svm)
    return obs_index


def get_svm_exact(all_svm, obs_svm, all_svm_ids, sub_mat_proba, norm="ranked", get_mut_rate=True):
    # Sort deltas by value
    # Sort the SVM values and keep the nuc_changes in the same order
    sorted_svm, sorted_changes = zip(*sorted(zip(all_svm, all_svm_ids)))

    mut_rates = []
    if get_mut_rate:
        for nuc_change in sorted_changes:
            pos, source, target = nuc_change.split(":")
            mut_rates.append(sub_mat_proba[source][target])
        assert not np.isnan(mut_rates).any()
        assert len(mut_rates) == len(sorted_svm)

    deltas_neg = [x for x in sorted_svm if x < 0]
    deltas_pos = [x for x in sorted_svm if x >= 0]

    # Transform deltas to phenotype: min=0, max=1, midpoint of the distribution = 0.5 (for delta = 0)
    if norm == "ranked":
        all_phenotype = [0. + 0.5 * (i + 1) / (len(deltas_neg) + 1) for i in range(len(deltas_neg))]
        all_phenotype += [0.5 + 0.5 * (i + 1) / (len(deltas_pos) + 1) for i in range(len(deltas_pos))]
        assert np.all(np.diff(all_phenotype) > 0)
    else:
        min_neg, max_neg = abs(min(deltas_neg)), 0.
        all_phenotype = [0. + 0.5 * (abs(x) - min_neg) / (max_neg - min_neg) for x in deltas_neg]
        min_pos, max_pos = 0., max(deltas_pos)
        all_phenotype += [0.5 + 0.5 * (x - min_pos) / (max_pos - min_pos)for x in deltas_pos]

        print("Before: First", all_phenotype[0], "Second:", all_phenotype[1], "Third:", all_phenotype[2])
        # Avoid 0 and 1 as Beta is defined on ]0, 1[
        all_phenotype = np.clip(all_phenotype, 1e-10, 1 - 1e-10)

    assert len(all_phenotype) == len(sorted_svm)
    assert min(all_phenotype) > 0.
    assert max(all_phenotype) < 1.

    obs_index = get_obs_index(obs_svm, sorted_svm)

    return np.array(mut_rates), obs_index, np.array(all_phenotype)

# Beta distribution probability density function
@nb.jit(nopython=True)
def beta_distribution(x: float, a: float, b: float) -> float:
    inv_beta = np.math.gamma(a + b) / (np.math.gamma(a) * np.math.gamma(b))
    return np.power(x, a - 1) * np.power(1 - x, b - 1) * inv_beta


# Selection coefficient from delta's quantile and model's parameters
@nb.jit(nopython=True)
def coeff_selection(bin_val: float, params: np.ndarray) -> float:
    assert 0 < bin_val < 1  # bin value needs to be between 0 and 1
    alpha = params[0] if len(params) > 0 else 1.0
    beta = params[1] if len(params) == 2 else alpha
    w_mutant = beta_distribution(bin_val, a=alpha, b=beta)
    w_ancestral = beta_distribution(0.5, a=alpha, b=beta)
    if w_mutant < 1.e-10 or w_ancestral < 1.e-50:
        return -np.infty
    s = np.log(w_mutant / w_ancestral)
    return s


# Scaled fixation probability from selection coefficient
@nb.jit(nopython=True)
def proba_fixation(s: float) -> float:
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


# Get probability of substitution for each bin: P(Mut) * P(Fix)
@nb.jit(nopython=True)
def proba_substitution(params: np.ndarray, mutations_proba: np.ndarray, scaled_bins: np.ndarray) -> np.ndarray:
    output_array = np.zeros(len(scaled_bins))
    for b in range(len(scaled_bins)):
        # Probability of mutation
        proba_mut = mutations_proba[b]

        # Probability of fixation
        s = coeff_selection(scaled_bins[b], params)
        proba_fix = proba_fixation(s)

        output_array[b] = proba_mut * proba_fix

    # Scaled output
    sum_output = np.sum(output_array)
    if sum_output == 0:
        return output_array
    else:
        return output_array / sum_output

@nb.jit(nopython=True)
def loglikelihood(mutations_proba: np.ndarray, obs_bins: np.ndarray, scaled_bins: np.ndarray, params: np.ndarray) -> float:
    if len(params) == 1 and params[0] < 1e-10:
        return -np.infty
    if len(params) == 2 and (params[0] < 1e-10 or params[1] < 1e-10):
        return -np.infty

    subs_proba = proba_substitution(params, mutations_proba, scaled_bins)
    models_lk = np.array([subs_proba[i] if i >= 0 else subs_proba[0] for i in obs_bins])

    if min(models_lk) <= 0.0:  # Avoid log(0)
        return -np.infty
    lnl = np.sum(np.log(models_lk))

    return lnl


def conclusion_purif(alpha):
    return "Stabilizing" if alpha > 1.0 else "Disruptive"


def conclusion_pos(alpha, beta):
    m = alpha/(alpha+beta)
    if m > 0.5:
        return f"Directional (+)"
    if m < 0.5:
        return f"Directional (-)"


def run_estimations(all_svm, all_svm_id, obs_svm, sub_mat_proba, alpha_threshold=0.05, min_bin=100, bins="hist", verbose=False):
    # Check if there are at least 2 different substitutions
    if len(set(obs_svm)) < 2:
        return None, None

    # Get the bin of the SVM distribution
    if bins == "quantile":
        mutations_proba, obs_bins, scaled_bins = get_svm_quantiles(all_svm, obs_svm, all_svm_id, sub_mat_proba, n_quant=min_bin)
    elif bins == "hist":
        mutations_proba, obs_bins, scaled_bins = get_svm_hist(all_svm, obs_svm, n_bin=min_bin)
    elif bins == "exact_ranked":
        mutations_proba, obs_bins, scaled_bins = get_svm_exact(all_svm, obs_svm, all_svm_id, sub_mat_proba, norm="ranked")
    elif bins == "exact_absolute":
        mutations_proba, obs_bins, scaled_bins = get_svm_exact(all_svm, obs_svm, all_svm_id, sub_mat_proba, norm="absolute")
    else:
        raise ValueError("Bins method should be in: 'quantile', 'hist', 'exact_ranked', 'exact_absolute'")

    total_bins = len(mutations_proba)

    # Null model: no param.
    ll_neutral = loglikelihood(mutations_proba, obs_bins, scaled_bins, np.array([]))

    # Stabilizing selection: alpha = beta
    bounds = [(0.0, np.inf)]
    model_purif = minimize(lambda theta: -loglikelihood(mutations_proba, obs_bins, scaled_bins, theta), np.array([1.0]),
                           bounds=bounds, method="Nelder-Mead") # L-BFGS-B "Nelder-Mead"
    ll_purif = -model_purif.fun

    # Positive selection
    bounds = [(0.0, np.inf), (0.0, np.inf)]
    initial_guess = np.array([1.0, 1.0])
    model_pos = minimize(lambda theta: -loglikelihood(mutations_proba, obs_bins, scaled_bins, theta), initial_guess,
                         bounds=bounds, method="Nelder-Mead")
    ll_pos = -model_pos.fun

    models = [model_purif, model_pos, bounds]

    # Likelihood ratio test:
    lrt_null_purif = -2 * (ll_neutral - ll_purif)
    lrt_purif_pos = -2 * (ll_purif - ll_pos)
    lrt_null_pos = -2 * (ll_neutral - ll_pos)

    p_value_null_purif = chi2.sf(lrt_null_purif, 1)
    p_value_purif_pos = chi2.sf(lrt_purif_pos, 1)
    p_value_null_pos = chi2.sf(lrt_null_pos, 2)

    conclusion = "Neutral model"
    # Purifying selection
    if p_value_null_purif < alpha_threshold:
        conclusion = conclusion_purif(model_purif.x[0])

    # Directional selection
    if p_value_null_pos < alpha_threshold and p_value_purif_pos < alpha_threshold:
        conclusion = conclusion_pos(model_pos.x[0], model_pos.x[1])

    # Low directional selection
    if p_value_null_pos < alpha_threshold < p_value_null_purif:
        conclusion = conclusion_pos(model_pos.x[0], model_pos.x[1])

    if verbose:
        print(f"  LnL neutral    : {ll_neutral:.2f}")
        print(f"  LnL stabilizing: {ll_purif:.4f}, alpha: {model_purif.x[0]:.3g}, beta: {model_purif.x[0]:.3g}")
        print(f"  LnL positive   : {ll_pos:.4f}, alpha: {model_pos.x[0]:.3g}, beta: {model_pos.x[1]:.3g}")
        print(f"  LRT null vs purif: {lrt_null_purif:.2f}, p-value: {p_value_null_purif:.2g}")
        print(f"  LRT purif vs pos : {lrt_purif_pos:.2f}, p-value: {p_value_purif_pos:.2g}")
        print("  Conclusion:", conclusion)

    # Create a DataFrame
    result = pd.DataFrame({"Nmut": [len(obs_svm)], "SumObs": [np.sum(obs_svm)], "MeanObs": [np.mean(obs_svm)],
                           "VarObs": [np.var(obs_svm)], "MedSVM": [np.median(all_svm)],
                           "MinSVM": [np.min(all_svm)], "MaxSVM": [np.max(all_svm)],
                           "AlphaPurif": [model_purif.x[0]], "AlphaPos": [model_pos.x[0]], "BetaPos": [model_pos.x[1]],
                           "NiterPurif": [model_purif.nit], "NiterPos": [model_pos.nit], "Nbins": [total_bins],
                           "LL_neutral": [ll_neutral], "LL_purif": [ll_purif], "LL_pos": [ll_pos],
                           "LRT_null_purif": [lrt_null_purif], "LRT_purif_pos": [lrt_purif_pos],
                           "p_value_null_purif": [p_value_null_purif], "p_value_purif_pos": [p_value_purif_pos],
                           "p_value_null_pos": [p_value_null_pos], "Conclusion": [conclusion]})

    return result, models


# Plot the results for each sequence
def general_plot(all_deltas, obs, svm_distribution, purif, pos, scenario, test, sub_ax1, sub_ax2):
    sub_ax1.hist(all_deltas, bins=100, density=True, alpha=0.7, label="All mutations")
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
    sub_ax1.set_title(f'{scenario}')

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
def plot_model(obs, all_svm, model_params, ax, model_type="Stabilizing", bins="hist", min_bin=50, bounds=None):
    if bins == "quantile":
        mutations_proba, obs_bins, scaled_bins = get_svm_quantiles(all_svm, obs, n_quant=min_bin)
    elif bins == "hist":
        mutations_proba, obs_bins, scaled_bins = get_svm_hist(all_svm, obs, n_bin=min_bin)

    k = 0.5
    if model_type == "Stabilizing":
        alpha_min = max([bounds[0][0], model_params[0] * k])
        alpha_max = min([bounds[0][1], model_params[0] * (2 - k)])
        alpha_range = np.linspace(alpha_min, alpha_max, 100)  # Std range
        ll_values = [loglikelihood(mutations_proba, obs_bins, scaled_bins, [std]) for std in alpha_range]
        ax.plot(alpha_range, ll_values, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], loglikelihood(mutations_proba, obs_bins, scaled_bins, model_params), color='red', marker='o',
                   label="Minimized Point")
        ax.set_xlabel("Parameter: Standard Deviation")
        ax.set_ylabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
        ax.legend()

    elif model_type == "Positive":
        alpha_min = max([bounds[0][0], model_params[0] * k])
        alpha_max = min([bounds[0][1], model_params[0] * (2 - k)])
        beta_min = max([bounds[1][0], model_params[1] * k])
        beta_max = min([bounds[1][1], model_params[1] * (2 - k)])
        alpha_range = np.linspace(alpha_min, alpha_max, 30)
        beta_range = np.linspace(beta_min, beta_max, 30)

        alpha_values, beta_values = np.meshgrid(alpha_range, beta_range)
        ll_values = np.array(
            [[loglikelihood(mutations_proba, obs_bins, scaled_bins, [std, mean]) for std in alpha_range] for mean in beta_range])

        # Plotting the log-likelihood surface
        ax.plot_surface(alpha_values, beta_values, ll_values, cmap='viridis', alpha=0.8,
                        label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], model_params[1], loglikelihood(mutations_proba, obs_bins, scaled_bins, model_params), color='red',
                   marker='o', label="Minimized Point")
        ax.invert_xaxis()
        ax.set_xlabel("Parameter: $\\alpha$")
        ax.set_ylabel("Parameter: $\\beta$")
        ax.set_zlabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")


'''
################################# !!!!!!! Temporary to test manually !!!!!!! ###########################################
DeltaSVM = pd.read_csv("ancestral_all_possible_deltaSVM.txt", sep='\t', header=0)
AllObsSVM = pd.read_csv("ancestral_to_observed_deltaSVM.txt", sep='\t', header=None, names=range(150+4))
AllObsSVM.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(AllObsSVM.columns[4:])

#ID = "chr13:86947569:86947706_13:86947569:86947706:Interval_6751"
result_list = []
for ID in AllObsSVM['ID']:
    all_svm_row = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0]
    obs_svm_row = AllObsSVM.loc[AllObsSVM['ID'] == ID, 4:].iloc[0]
    all_svm = all_svm_row.dropna().values.tolist()
    obs_svm = obs_svm_row.dropna().values.tolist()
    results, model = run_estimations(all_svm, obs_svm, alpha_threshold=0.01, min_bin=50, bins="quantile", verbose=False)

    result_list.append(results)

#final_results = pd.concat(result_list)
#final_results.to_csv(f"MLE_results_test_quantile.csv", index=False)
########################################################################################################################
'''
