#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import pandas as pd
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from alive_progress import alive_bar
from multiprocessing import Pool
import numba as nb

GREEN = "#8FB03E"
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
exp_name = "lognorm_gibbs_nb"

np.random.seed(0)


# Beta distribution probability density function
@nb.jit(nb.float64(nb.float64, nb.float64, nb.float64), nopython=True)
def beta_distribution(x, a, b):
    beta_val = np.math.gamma(a + b) / (np.math.gamma(a) * np.math.gamma(b))
    return np.power(x, a - 1) * np.power(1 - x, b - 1) / beta_val


# Selection coefficient from delta's quantile and model's parameters
@nb.jit(nb.float64(nb.float64, nb.float64, nb.float64), nopython=True)
def coeff_selection(quant_delta, alpha, beta):
    assert 0 < quant_delta < 1
    w_mutant = beta_distribution(quant_delta, a=alpha, b=beta)
    w_ancestral = beta_distribution(0.5, a=alpha, b=beta)
    if w_mutant < 1.e-10 or w_ancestral < 1.e-10:
        return -np.infty
    s = np.log(w_mutant / w_ancestral)
    return s


# Scaled fixation probability from selection coefficient
@nb.jit(nb.float64(nb.float64), nopython=True)
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
@nb.jit(nb.float64[:](nb.float64, nb.float64, nb.float64[:]), nopython=True)
def proba_substitution(alpha, beta, mutations_proba):
    n_bins = len(mutations_proba)
    output_array = np.zeros(n_bins)
    delta = 1 / (2 * n_bins)
    for quant in range(n_bins):
        # Probability of mutation
        proba_mut = mutations_proba[quant]

        # Probability of fixation
        quant_val = quant / n_bins + delta  # quantile needs to be between 0<quant<1
        s = coeff_selection(quant_val, alpha, beta)
        proba_fix = proba_fixation(s)

        output_array[quant] = proba_mut * proba_fix

    # Scaled output
    sum_output = np.sum(output_array)
    if sum_output == 0:
        return output_array
    else:
        return output_array / sum_output


# Log likelihood of the model
@nb.jit(nb.float64(nb.float64[:], nb.int64[:], nb.float64, nb.float64), nopython=True)
def loglikelihood(mutations_proba, obs_bins, alpha, beta):
    subs_proba = proba_substitution(alpha, beta, mutations_proba)
    models_lk = np.array([subs_proba[i] if i >= 0 else subs_proba[0] for i in obs_bins])

    if np.min(models_lk) <= 0.0:  # Avoid log(0)
        return -np.infty
    lnl = np.sum(np.log(models_lk))
    return lnl


# Normal distribution in log scale
@nb.jit(nb.float64(nb.float64), nopython=True)
def norm_distribution_log(x):
    # normal distribution in log scale
    return -0.5 * np.log(2 * np.pi) - 0.5 * np.power(x, 2)


# Log scaling parameter
@nb.jit(nb.float64(nb.float64), nopython=True)
def log_scaling_param(x):
    # normal distribution in log scale, with x = log10(x)
    return norm_distribution_log(np.log10(x))


# Log prior of the model
@nb.jit(nb.float64(nb.float64, nb.float64), nopython=True)
def logprior(alpha, beta):
    return log_scaling_param(alpha) + log_scaling_param(beta)


# Move alpha parameter with Metropolis-Hastings
@nb.jit(nb.types.UniTuple(nb.float64, 3)(
    nb.float64, nb.float64, nb.float64[:], nb.int64[:], nb.float64, nb.float64, nb.float64), nopython=True)
def move_alpha(lnp, lnl, mutations_proba, obs_bins, alpha, beta, tuning):
    m = tuning * (np.random.uniform() - 0.5)
    new_alpha = np.exp(m) * alpha
    new_lnp = logprior(new_alpha, beta)
    new_lnl = loglikelihood(mutations_proba, obs_bins, new_alpha, beta)
    delta_logprob = (new_lnp + new_lnl) - (lnp + lnl) + m
    accepted = (np.log(np.random.uniform()) < delta_logprob)
    if accepted:
        return new_alpha, new_lnp, new_lnl
    else:
        return alpha, lnp, lnl


# Move beta parameter with Metropolis-Hastings
@nb.jit(nb.types.UniTuple(nb.float64, 3)(
    nb.float64, nb.float64, nb.float64[:], nb.int64[:], nb.float64, nb.float64, nb.float64), nopython=True)
def move_beta(lnp, lnl, mutations_proba, obs_bins, alpha, beta, tuning):
    m = tuning * (np.random.uniform() - 0.5)
    new_beta = np.exp(m) * beta
    new_lnp = logprior(alpha, new_beta)
    new_lnl = loglikelihood(mutations_proba, obs_bins, alpha, new_beta)
    delta_logprob = (new_lnp + new_lnl) - (lnp + lnl) + m
    accepted = (np.log(np.random.uniform()) < delta_logprob)
    if accepted:
        return new_beta, new_lnp, new_lnl
    else:
        return beta, lnp, lnl


@nb.jit(nb.types.UniTuple(nb.float64, 4)(
    nb.float64[:], nb.int64[:], nb.float64, nb.float64, nb.float64, nb.float64, nb.float64), nopython=True)
def move_gibbs(mutations_proba, obs_bins, alpha, beta, lnp, lnl, tuning):
    iter_per_point = 10
    for i in range(iter_per_point):
        alpha, lnp, lnl = move_alpha(lnp, lnl, mutations_proba, obs_bins, alpha, beta, tuning)
        beta, lnp, lnl = move_beta(lnp, lnl, mutations_proba, obs_bins, alpha, beta, tuning)
    return alpha, beta, lnp, lnl


@nb.jit(nb.types.UniTuple(nb.float64[:, :], 2)(
    nb.float64[:], nb.int64[:], nb.int64, nb.float64, nb.float64), nopython=True)
def mcmc(mutations_proba, obs_bins, n_points, alpha, beta):
    lnp = logprior(alpha, beta)
    lnl = loglikelihood(mutations_proba, obs_bins, alpha, beta)
    # Initialize the chain
    chain_params = [[alpha, beta]]
    chain_logprob = [[lnp, lnl]]

    # Run the MCMC
    for i in range(n_points):
        alpha, beta, lnp, lnl = move_gibbs(mutations_proba, obs_bins, alpha, beta, lnp, lnl, 0.05)
        alpha, beta, lnp, lnl = move_gibbs(mutations_proba, obs_bins, alpha, beta, lnp, lnl, 0.15)

        chain_params.append([alpha, beta])
        chain_logprob.append([lnp, lnl])

    chain_params = np.array(chain_params)
    chain_logprob = np.array(chain_logprob)
    return chain_params, chain_logprob


# Get the quantiles of the SVM distribution of each side
def get_svm_quantiles(all_svm, obs_svm, n_quant):
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
    return np.array(quant_proba), obs_bins, np.array(obs_quant)


def run_mcmc(all_svm, obs_svm, n_quant, n_points, alpha, beta):
    mutations_proba, obs_bins, obs_quant = get_svm_quantiles(all_svm, obs_svm, n_quant)
    chain_params, chain_logprob = mcmc(mutations_proba, obs_bins, n_points, alpha, beta)
    cdf_50 = stats.beta.cdf(0.5, a=chain_params[:, 0], b=chain_params[:, 1])
    chain_params = np.column_stack((chain_params, cdf_50))
    return chain_params, chain_logprob


def conclusion(alpha, beta):
    if alpha > 1.0 and beta > 1.0:
        return "Stabilizing " + ("(+)" if alpha > beta else "(-)")
    elif alpha > 1.0 > beta:
        return "Directional (+)"
    elif alpha < 1.0 < beta:
        return "Directional (-)"
    elif alpha < 1.0 and beta < 1.0:
        return "Disruptive " + ("(+)" if alpha > beta else "(-)")
    else:
        assert alpha == 1.0 and beta == 1.0
        return "Neutral"


def chain_results(input_chain_params, burnin):
    # Remove burnin
    n_burnin = int(burnin * len(input_chain_params))
    chain_params = input_chain_params[n_burnin:]
    post_mean_alpha, post_mean_beta, post_mean_cdf_50 = np.mean(chain_params, axis=0)
    # print(f"Posterior mean alpha: {post_mean_alpha:.2f}, beta: {post_mean_beta:.2f}")
    chain_beta_exp = chain_params[:, 0] / (chain_params[:, 0] + chain_params[:, 1])
    post_mean_beta_exp = np.mean(chain_beta_exp)
    # print(f"Posterior mean beta_exp: {post_mean_beta_exp:.2f}")
    # print(f"Posterior mean cdf_50: {post_mean_cdf_50:.2f}")
    conclusion_list = [conclusion(alpha, beta) for alpha, beta, cdf_50 in chain_params]
    conclusion_count = Counter(conclusion_list)
    conclusion_proba = {k: conclusion_count[k] / len(conclusion_list) for k in sorted(conclusion_count.keys())}
    # print(f"Conclusions: {conclusion_proba}")
    return post_mean_alpha, post_mean_beta, post_mean_beta_exp, post_mean_cdf_50, conclusion_proba


def save_trace(chain_params, chain_logprob, ID):
    trace_dict = defaultdict(list)
    trace_dict["alpha"].extend(chain_params[:, 0])
    trace_dict["beta"].extend(chain_params[:, 1])
    trace_dict["logprior"].extend(chain_logprob[:, 0])
    trace_dict["loglikelihood"].extend(chain_logprob[:, 1])
    trace_dict["logposterior"].extend(chain_logprob[:, 0] + chain_logprob[:, 1])
    trace_df = pd.DataFrame(trace_dict)
    trace_df.to_csv(f"trace_Bayes_{exp_name}/{ID}_trace.csv.gz", index=False)


def plot_trace(chain_params_list, chain_logprob_list, ID, burnin):
    n_burnin = int(burnin * len(chain_params_list[0]))
    fig, axes = plt.subplots(2, 2, figsize=(16, 9))
    colors = [GREEN, RED, YELLOW, BLUE, LIGHTGREEN]
    for i, chain_params in enumerate(chain_params_list):
        color = colors[i % len(colors)]
        axes[0, 0].plot(chain_params[:, 0], label=f"Chain {i + 1} - α", color=color, linestyle="--")
        axes[0, 0].plot(chain_params[:, 1], label=f"Chain {i + 1} - β", color=color, linestyle="-")
    axes[0, 0].axhline(y=1.0, color='black', linestyle='--')
    axes[0, 0].axvline(x=n_burnin, color='black', linestyle='-')
    axes[0, 0].set_title("Trace of α and β")
    axes[0, 0].legend()
    for i, chain_logprob in enumerate(chain_logprob_list):
        color = colors[i % len(colors)]
        axes[0, 1].plot(chain_logprob[:, 0], label=f"Chain {i + 1}", color=color)
        axes[1, 0].plot(chain_logprob[:, 1], label=f"Chain {i + 1}", color=color)
        axes[1, 1].plot(chain_logprob[:, 0] + chain_logprob[:, 1], label=f"Chain {i + 1}", color=color)
    axes[0, 1].set_title("Log prior")
    axes[1, 0].set_title("Log likelihood")
    axes[1, 1].set_title("Log posterior")
    axes[0, 1].axvline(x=n_burnin, color='black', linestyle='-')
    axes[1, 0].axvline(x=n_burnin, color='black', linestyle='-')
    axes[1, 1].axvline(x=n_burnin, color='black', linestyle='-')
    axes[0, 1].legend()
    axes[1, 0].legend()
    axes[1, 1].legend()
    plt.savefig(f"plots_Bayes_{exp_name}/{ID}_trace.pdf")
    plt.close()


def run_estimations(all_svm, obs_svm, ID, alpha, beta, burnin):
    n_points = 10000
    n_quant = 25
    chain_params, chain_logprob = run_mcmc(all_svm, obs_svm, n_quant, n_points, alpha, beta)
    save_trace(chain_params, chain_logprob, ID)
    pm_alpha, pm_beta, pm_beta_exp, pm_cdf_50, conclusion_proba = chain_results(chain_params, burnin)
    dico_results = {"ID": [ID], "alpha": [pm_alpha], "beta": [pm_beta], "m": [pm_beta_exp], "cdf_50": [pm_cdf_50]}
    for k, v in conclusion_proba.items():
        dico_results[k] = [v]
    return pd.DataFrame(dico_results), chain_logprob, chain_params


################################# !!!!!!! Temporary to test manually !!!!!!! ###########################################
DeltaSVM = pd.read_csv("ancestral_all_possible_deltaSVM.txt", sep='\t', header=0)
AllObsSVM = pd.read_csv("ancestral_to_observed_deltaSVM.txt", sep='\t', header=None, names=range(150 + 4))
AllObsSVM.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(AllObsSVM.columns[4:])


def get_row_estimate(ID, burnin=0.5):
    all_svm_row = DeltaSVM.loc[DeltaSVM['ID'] == ID, "pos0:A":].iloc[0].dropna().values.tolist()
    obs_svm_row = AllObsSVM.loc[AllObsSVM['ID'] == ID, 4:].iloc[0].dropna().values.tolist()
    # print(f"Obs SVM:{obs_svm_row}")
    results_1, cl_1, cp_1 = run_estimations(all_svm_row, obs_svm_row, f"{ID}_1", 0.1, 0.1, burnin)
    results_2, cl_2, cp_2 = run_estimations(all_svm_row, obs_svm_row, f"{ID}_2", 4.0, 4.0, burnin)
    results_3, cl_3, cp_3 = run_estimations(all_svm_row, obs_svm_row, f"{ID}_3", 1.0, 1.0, burnin)
    plot_trace([cp_1, cp_2, cp_3], [cl_1, cl_2, cl_3], ID, burnin)
    return results_1, results_2, results_3


os.makedirs(f"plots_Bayes_{exp_name}", exist_ok=True)
os.makedirs(f"trace_Bayes_{exp_name}", exist_ok=True)
NbThread = 4
NRows = 100
# ID = "chr13:86947569:86947706_13:86947569:86947706:Interval_6751"
if __name__ == '__main__':
    result_list = []
    with alive_bar(NRows) as bar:  # progress bar
        with Pool(NbThread) as pool:
            # Run function for each sequence in parallel
            for results in pool.imap_unordered(get_row_estimate, AllObsSVM['ID'][:NRows]):
                for res in results:
                    result_list.append(res)
                bar()
    final_results = pd.concat(result_list, axis=0)
    final_results.fillna(0.0, inplace=True)
    final_results.sort_values(by='ID', inplace=True)
    final_results.to_csv(f"Bayesian_results_{exp_name}.csv", index=False)
########################################################################################################################
