#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


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
    return quant_proba, obs_bins, obs_quant


# Selection coefficient from delta's quantile and model's parameters
def coeff_selection(quant_delta, alpha, beta):
    assert 0 < quant_delta < 1
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


def loglikelihood(mutations_proba, obs_bins, alpha, beta):
    subs_proba = proba_substitution(alpha, beta, mutations_proba)
    models_lk = [subs_proba[i] if i >= 0 else subs_proba[0] for i in obs_bins]

    if min(models_lk) <= 0.0:  # Avoid log(0)
        return -np.infty
    lnl = np.sum(np.log(models_lk))
    return lnl


def logprior(alpha, beta):
    # Add a penalty params far from 1.0
    # lnp = stats.gamma.logpdf(alpha, a=2, scale=0.5) + stats.gamma.logpdf(beta, a=2, scale=0.5)
    lnp = abs(np.log(alpha)) + abs(np.log(beta))
    return -lnp


def logprob(mutations_proba, obs_bins, alpha, beta):
    return logprior(alpha, beta), loglikelihood(mutations_proba, obs_bins, alpha, beta),


def move_param(logprob_value, logprob_func, param, tuning):
    m = tuning * (np.random.uniform() - 0.5)
    new_param = np.exp(m) * param
    new_logprob = logprob_func(new_param)

    deltalogprob = sum(new_logprob) - sum(logprob_value) + m
    accepted = (np.log(np.random.uniform()) < deltalogprob)
    if accepted:
        return new_param, new_logprob
    else:
        return param, logprob_value


def move_mcmc(mutations_proba, obs_bins, alpha, beta, logprob_value, tuning, iter_per_point):
    for i in range(iter_per_point):
        alpha, logprob_value = move_param(logprob_value,
                                          lambda x: logprob(mutations_proba, obs_bins, x, beta), alpha, tuning)
        beta, logprob_value = move_param(logprob_value,
                                         lambda x: logprob(mutations_proba, obs_bins, alpha, x), beta, tuning)
    return alpha, beta, logprob_value


def run_mcmc(all_svm, obs_svm, n_quant, n_points, iter_per_point, tuning, alpha, beta):
    mutations_proba, obs_bins, obs_quant = get_svm_quantiles(all_svm, obs_svm, n_quant)

    logprob_value = logprob(mutations_proba, obs_bins, alpha, beta)
    # Initialize the chain
    chain_params = [(alpha, beta)]
    chain_logprob = [logprob_value]

    # Run the MCMC
    for i in range(n_points):
        alpha, beta, logprob_value = move_mcmc(mutations_proba, obs_bins, alpha, beta, logprob_value, tuning,
                                               iter_per_point)
        chain_params.append((alpha, beta))
        chain_logprob.append(logprob_value)

    return np.array(chain_params), np.array(chain_logprob)


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
    post_mean_alpha, post_mean_beta = np.mean(chain_params, axis=0)
    print(f"Posterior mean alpha: {post_mean_alpha:.2f}, beta: {post_mean_beta:.2f}")
    chain_beta_exp = chain_params[:, 0] / (chain_params[:, 0] + chain_params[:, 1])
    post_mean_beta_exp = np.mean(chain_beta_exp)
    print(f"Posterior mean beta_exp: {post_mean_beta_exp:.2f}")
    conclusion_list = [conclusion(alpha, beta) for alpha, beta in chain_params]
    conclusion_count = Counter(conclusion_list)
    conclusion_proba = {k: conclusion_count[k] / len(conclusion_list) for k in sorted(conclusion_count.keys())}
    print(f"Conclusions: {conclusion_proba}")
    return post_mean_alpha, post_mean_beta, post_mean_beta_exp, conclusion_proba


def plot_trace(chain_params, chain_logprob, ID, burnin):
    n_burnin = int(burnin * len(chain_params))
    fig, axes = plt.subplots(2, 2, figsize=(16, 9))
    axes[0, 0].plot(chain_params[:, 0], label="alpha")
    axes[0, 0].plot(chain_params[:, 1], label="beta")
    axes[0, 0].axhline(y=1.0, color='r', linestyle='--')
    axes[0, 0].axvline(x=n_burnin, color='black', linestyle='-')
    axes[0, 0].set_title("Trace of alpha and beta")
    axes[0, 0].legend()
    axes[0, 1].plot(chain_logprob[:, 0])
    axes[0, 1].axvline(x=n_burnin, color='black', linestyle='-')
    axes[0, 1].set_title("Log prior")
    axes[1, 0].plot(chain_logprob[:, 1])
    axes[1, 0].axvline(x=n_burnin, color='black', linestyle='-')
    axes[1, 0].set_title("Log likelihood")
    axes[1, 1].plot(chain_logprob[:, 0] + chain_logprob[:, 1])
    axes[1, 1].axvline(x=n_burnin, color='black', linestyle='-')
    axes[1, 1].set_title("Log posterior")
    plt.savefig(f"plots_Bayes/{ID}_trace.png")
    plt.close()


def run_estimations(all_svm, obs_svm, ID, alpha, beta):
    n_points, iter_per_point = 20000, 1
    n_quant = 25
    tuning = 0.1
    burnin = 0.5
    chain_params, chain_logprob = run_mcmc(all_svm, obs_svm, n_quant, n_points, iter_per_point, tuning, alpha, beta)
    plot_trace(chain_params, chain_logprob, ID, burnin)
    post_mean_alpha, post_mean_beta, post_mean_beta_exp, conclusion_proba = chain_results(chain_params, burnin)
    dico_results = {"ID": [ID], "alpha": [post_mean_alpha], "beta": [post_mean_beta], "m": [post_mean_beta_exp]}
    for k, v in conclusion_proba.items():
        dico_results[k] = [v]
    return pd.DataFrame(dico_results)


################################# !!!!!!! Temporary to test manually !!!!!!! ###########################################
DeltaSVM = pd.read_csv("ancestral_all_possible_deltaSVM.txt", sep='\t', header=0)
AllObsSVM = pd.read_csv("ancestral_to_observed_deltaSVM.txt", sep='\t', header=None, names=range(150 + 4))
AllObsSVM.columns = ['ID', 'SVM', 'Total_deltaSVM', 'NbSub'] + list(AllObsSVM.columns[4:])

# ID = "chr13:86947569:86947706_13:86947569:86947706:Interval_6751"
os.makedirs(f"plots_Bayes", exist_ok=True)
result_list = []
for ID_row in AllObsSVM['ID'][:20]:
    all_svm_row = DeltaSVM.loc[DeltaSVM['ID'] == ID_row, "pos0:A":].iloc[0].dropna().values.tolist()
    obs_svm_row = AllObsSVM.loc[AllObsSVM['ID'] == ID_row, 4:].iloc[0].dropna().values.tolist()
    print(f"Obs SVM:{obs_svm_row}")
    results_1 = run_estimations(all_svm_row, obs_svm_row, f"{ID_row}_1", 0.1, 0.1)
    result_list.append(results_1)
    results_2 = run_estimations(all_svm_row, obs_svm_row, f"{ID_row}_2", 10.0, 10.0)
    result_list.append(results_2)

final_results = pd.concat(result_list, axis=0)
final_results.fillna(0.0, inplace=True)
final_results.to_csv(f"Bayesian_results.csv", index=False)
########################################################################################################################
