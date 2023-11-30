#!/usr/bin/env python
# coding=utf-8
import numpy as np
from scipy import stats
import matplotlib as mat


def proba_fixation(delta, params, delta_bounds):
    alpha = params[0] if len(params) > 0 else 1.0
    beta = params[1] if len(params) == 2 else alpha
    max_delta = max(np.abs(delta_bounds))
    scaled_delta = (delta + max_delta) / (2 * max_delta)
    w_mutant = stats.beta.pdf(scaled_delta, a=alpha, b=beta)
    w_ancestral = stats.beta.pdf(0.5, a=alpha, b=beta)
    if w_mutant < 1.e-10:
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
    delta_bounds = [min(bins_values), max(bins_values)]
    output_array = np.zeros(n_bins)
    for b in range(n_bins):
        p = mutations_proba[b]
        min_bin = bins_values[b]
        max_bin = bins_values[b + 1]
        mean_bin = (max_bin + min_bin) / 2
        output_array[b] = p * proba_fixation(mean_bin, params, delta_bounds)
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
def plot_model(obs, hist_svm, model_params, ax, model_type="Stabilizing", bounds=None):
    k = 0.5
    if model_type == "Stabilizing":
        alpha_min = max([bounds[0][0], model_params[0] * k])
        alpha_max = min([bounds[0][1], model_params[0] * (2 - k)])
        alpha_range = np.linspace(alpha_min, alpha_max, 100)  # Std range
        ll_values = [-loglikelihood(obs, [std], hist_svm) for std in alpha_range]
        ax.plot(alpha_range, ll_values, label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], -loglikelihood(obs, model_params, hist_svm), color='red', marker='o',
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
            [[-loglikelihood(obs, [std, mean], hist_svm) for std in alpha_range] for mean in beta_range])

        # Plotting the log-likelihood surface
        ax.plot_surface(alpha_values, beta_values, ll_values, cmap='viridis', alpha=0.8,
                        label=f"{model_type} Selection Model")
        ax.scatter(model_params[0], model_params[1], -loglikelihood(obs, model_params, hist_svm), color='red',
                   marker='o', label="Minimized Point")
        ax.invert_xaxis()
        ax.set_xlabel("Parameter: $\\alpha$")
        ax.set_ylabel("Parameter: $\\beta$")
        ax.set_zlabel("Log-Likelihood")
        ax.set_title(f"Minimization Process - {model_type} Selection Model")
