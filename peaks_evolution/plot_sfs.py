import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

sfs_weight = {"watterson": lambda i, n: 1.0 / i, "tajima": lambda i, n: n - i, "fay_wu": lambda i, n: i}

fontsize = 16
fontsize_legend = 14
my_dpi = 128
GREEN = "#8FB03E"
RED = "#EB6231"
YELLOW = "#E29D26"
BLUE = "#5D80B4"
LIGHTGREEN = "#6ABD9B"
LIGHTYELLOW = "#FFFCBF"


def theta(sfs_epsilon, daf_n, weight_method):
    sfs_theta = sfs_epsilon * np.array(range(1, daf_n))
    weights = np.array([sfs_weight[weight_method](i, daf_n) for i in range(1, daf_n)])
    return sum(sfs_theta * weights) / sum(weights)


def daf_to_sfs(daf_list, min_n):
    array_daf = np.array(daf_list, dtype=np.int64).T
    return np.array([np.bincount(d, minlength=min_n) for d in array_daf])


def normalize_sfs(sfs):
    return (sfs.T / np.sum(sfs, axis=1)).T


def plot_sfs(snp_sfs, max_daf, daf_axis, scaled, cat_mean, cat_color, ax=None):
    min_mean, max_mean = np.inf, -np.inf
    for cat in snp_sfs:
        sfs = snp_sfs[cat][:, 1:].copy()
        if scaled == "neutral":
            sfs *= np.array([i for i in range(1, max_daf)])
        elif scaled == "normalize":
            sfs = normalize_sfs(sfs)

        mean_sfs = np.mean(sfs, axis=0)
        min_mean = min(min_mean, np.min(mean_sfs))
        max_mean = max(max_mean, np.max(mean_sfs))
        std_sfs = np.std(sfs, axis=0)
        f_sfs = ((mean_sfs > 0) & np.isfinite(mean_sfs))
        ax.scatter(daf_axis[f_sfs], mean_sfs[f_sfs], s=1.0, color=cat_color[cat])
        ax.plot(daf_axis[f_sfs], mean_sfs[f_sfs], label=cat_mean[cat], linewidth=1.0, color=cat_color[cat])
        if len(sfs) > 1:
            ax.fill_between(daf_axis, mean_sfs - std_sfs, mean_sfs + std_sfs, linewidth=1.0, alpha=0.2,
                            color=cat_color[cat])
    if max_daf < 32:
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_xlabel("Derived allele count")
    if scaled != "neutral":
        ax.set_ylim(min_mean * 0.5, max_mean * 1.2)
        ax.set_yscale("log")
    ax.set_ylabel("Proportion of mutations")
    ax.set_title(f"SFS ({scaled})")
    ax.legend()


def AIC_weight(row):
    aic_neutral = -2 * row["LL_neutral"]
    aic_purif = -2 * row["LL_purif"] + 1 * 2
    aic_pos = -2 * row["LL_pos"] + 2 * 2
    max_aic = min(aic_neutral, aic_purif, aic_pos)
    weights = np.exp(-0.5 * (np.array([aic_neutral, aic_purif, aic_pos]) - max_aic))
    return weights / np.sum(weights)


def selcoef_compound(row):
    s_array = np.array([0, row["SelCoefStab"], row["SelCoefPos"]])
    w_array = np.array([row["w_neutral"], row["w_pur"], row["w_pos"]])
    filter_array = np.isfinite(s_array)
    # filter_array[-1] = False
    return np.sum(s_array[filter_array] * w_array[filter_array]) / np.sum(w_array[filter_array])


def transform_s(x):
    return x if x != 0.0 else np.nan


def main(output_pdf, input_mle, input_snp, subsample=24, nbr_replicates=100):
    output_dir = os.path.dirname(output_pdf)
    if not os.path.exists(output_dir) and output_dir != "":
        os.makedirs(os.path.dirname(output_pdf), exist_ok=True)

    df_mle = pd.read_csv(input_mle)
    df_mle[["w_neutral", "w_pur", "w_pos"]] = df_mle.apply(AIC_weight, axis=1, result_type="expand")

    df_snps = pd.read_csv(input_snp, sep='\t')
    df_snps = df_snps[df_snps["Flag"] == "ref_ancestral"]
    df_snps["SelCoefStab"] = df_snps["SelCoefStab"].apply(transform_s)
    df_snps["SelCoefPos"] = df_snps["SelCoefPos"].apply(transform_s)
    nn_nan = len(df_snps) - np.sum(np.isfinite(df_snps["SelCoefPos"]))
    print(f"Number of NaN in SelCoefPos: {nn_nan} / {len(df_snps)} = {nn_nan / len(df_snps) * 100:.2f}%")
    nn_nan = len(df_snps) - np.sum(np.isfinite(df_snps["SelCoefStab"]))
    print(f"Number of NaN in SelCoefStab: {nn_nan} / {len(df_snps)} = {nn_nan / len(df_snps) * 100:.2f}%")

    qcut = 5
    df = pd.merge(df_snps, df_mle, on="ID", how="left")
    assert len(df) == len(df_snps)
    assert len(set(df["NbTot"])) == 1
    df["SelCoefMixed"] = df.apply(selcoef_compound, axis=1)
    df = df[df["AlphaPurif"] > 1.0]
    df["cat"] = pd.qcut(df["SelCoefMixed"], qcut, duplicates="drop")

    # Use a red to blue color map
    cmap = plt.get_cmap("coolwarm", qcut)
    snps_daf = defaultdict(list)
    cat_mean, cat_color = {}, {}

    sample_size = max(df["NbTot"])
    max_daf = min(sample_size, subsample)
    gp = df.groupby("cat", observed=True)
    for i, (cat, df_cat) in enumerate(gp):
        # mask of finite values
        cat_mean[cat] = f"S={df_cat['SelCoefMixed'].mean():.3f} ({len(df_cat)} SNPs)"
        cat_color[cat] = cmap(i)
        for k in df_cat["NbAlt"]:
            if max_daf < sample_size:
                count = [np.random.hypergeometric(k, sample_size - k, max_daf) for _ in range(nbr_replicates)]
                count = [i if i != max_daf else 0 for i in count]
            else:
                count = [k]
            snps_daf[cat].append(count)

    snp_sfs = {cat: daf_to_sfs(daf, max_daf) for cat, daf in snps_daf.items()}
    snp_sfs_mean = {cat: np.mean(sfs, axis=0) for cat, sfs in snp_sfs.items()}

    theta_dict = defaultdict(list)
    daf_axis = np.array(range(1, max_daf))
    for cat, mean_sfs in snp_sfs_mean.items():
        theta_dict["category"].append(cat)
        theta_dict["cat_mean"].append(cat_mean[cat])
        theta_dict["daf_mean"].append(mean_sfs[1:].mean())
        for theta_method in sfs_weight:
            theta_dict[theta_method].append(theta(mean_sfs[1:], max_daf, theta_method))
    df_out = pd.DataFrame(theta_dict)
    print(df_out)
    df_out.to_csv(output_pdf.replace('.pdf', '.theta.csv'), index=False)

    fig, axes = plt.subplots(figsize=(1920 / my_dpi, 960 / my_dpi), dpi=my_dpi, nrows=2, ncols=2)
    for i, scaled in enumerate(["raw", "neutral", "normalize"]):
        plot_sfs(snp_sfs, max_daf, daf_axis, scaled, cat_mean, cat_color, axes[i // 2, i % 2])
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf")
    plt.clf()
    plt.close("all")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_MLE', required=False, type=str, dest="tsv", help="Input MLE file",
                        default="human_HNF6_MLE_summary.csv")
    parser.add_argument('--input_SNP', required=False, type=str, dest="snp", help="Input SNP file",
                        default="human_HNF6_SNP_SelectionCoefficient.txt")
    parser.add_argument('--output_pdf', required=False, type=str, dest="output_pdf", help="Output pdf file",
                        default="human_HNF6_SFS.pdf")
    args = parser.parse_args()
    main(args.output_pdf, args.tsv, args.snp)
