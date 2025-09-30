# Implement rules to test for Positive selection and simulate sequences evolution:
from snakemake.io import touch
import os

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
AncNode = config["AncNode"]
baseDir = os.path.abspath(config["baseDir"])

pathResults = f"{baseDir}/results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"{baseDir}/results/peaks_calling/{peakType}/{sp}/{sample}"

rule ComputeAllDeltaSVM:
    """Compute all possible and observed SVM"""
    wildcard_constraints: AncNode="(?!focal$)[a-zA-Z0-9_]+" # Exclude focal node
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_{AncNode}_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_{AncNode}_sequences_upper.fa"
    output:
        AllSVM = pathResults + "/{TF}/deltas/{AncNode}_all_possible_deltaSVM.txt",
        ObsSVM = pathResults + "/{TF}/deltas/{AncNode}_to_observed_deltaSVM.txt"
    log: out = pathResults + "/log/{TF}/ComputeAllDeltaSVM_{AncNode}.out"
    conda: "../../envs/selection_tests.yaml"
    threads: config["nbPart"]
    params: time="1:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python ../scripts/RegEvol/compute_all_deltaSVM.py {sp} {sample}/{wildcards.TF} {peakType} --node {AncNode} -T {threads} > {log.out} 2>&1 || exit 1
        """

def evaluate_time(ancestral_seq):
    # Calculate time needed to perform Permutations, considering 1500 peaks per hour and per thread.
    with open(ancestral_seq, 'r') as file:
        num_lines = sum(1 for _ in file)
        time = num_lines/ (1500*config["nbPart"])+2
        cluster_time = f"{round(time)}:00:00"
    return cluster_time

rule PermutationTest:
    message: "Test for positive selection between ancestral and focal sequences"
    input:
        AllSVM = pathResults + "/{TF}/deltas/{AncNode}_all_possible_deltaSVM.txt",
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_{AncNode}_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_{AncNode}_sequences_upper.fa"
    output: touch(pathResults + "/{TF}/Tests/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations_two_tailed_{AncNode}.txt")
    threads: max(10, config["nbPart"])
    log: out=pathResults + "/log/{TF}/PermutationTest_{AncNode}.out"
    params: time=lambda wildcards, input: evaluate_time(input.AllSVM), mem="5G", threads=max(10, config["nbPart"])
    conda: "../../envs/selection_tests.yaml"
    shell:
        """
        python ../scripts/permutations_test.py {sp} {sample} {wildcards.TF} --peakType {peakType} \
        --NbRand {config[nbRand]} --node {AncNode} --NbThread {threads} > {log.out} 2>&1 || exit 1
        """

rule MaxLLTest:
    message: "Test for positive selection using Likelihood Ratio Test between 3 Maximized models"
    input:
        AllSVM = pathResults + "/{TF}/deltas/{AncNode}_all_possible_deltaSVM.txt",
        ObsSVM = pathResults + "/{TF}/deltas/{AncNode}_to_observed_deltaSVM.txt"
    output: touch(pathResults + "/{TF}/Tests/MLE_summary_{BinType}_{AncNode}.csv")
    log: out = pathResults + "/log/{TF}/MaxLLTest_{BinType}_{AncNode}.out"
    threads: config["nbPart"]
    params: mem="10G", threads=config["nbPart"], nbBin=config["nbBin"], BinType=config["BinType"],
            threshold=config["threshold"], time="1:00:00"
    priority: 10
    conda: "../../envs/selection_tests.yaml"
    shell:
        """
        python ../scripts/RegEvol/MaxLL_estimation.py {sp} {sample}/{wildcards.TF} \
        --peakType {peakType} --binType {params.BinType} --NbBin {params.nbBin} --threshold {params.threshold} \
        --node {AncNode} -T {threads} > {log.out} 2>&1 
        """
