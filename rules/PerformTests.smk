# Implement rules to test for Positive selection and simulate sequences evolution:
from snakemake.io import touch

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"

rule PermutationTest:
    message: "Test for positive selection between ancestral and focal sequences"
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa"
    output: touch(pathResults + "/{TF}/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations.txt")
    threads: config["nbPart"]
    priority: 2
    log: out=pathResults + "/log/{TF}/PermutationTest.out"
    params: time="15:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python Positive_Selection_Tests/Permutation_Test/permutations.py {sp} {sample} {wildcards.TF} {peakType} \
        --NbRand {config[nbRand]} --{cluster} --NbThread {threads} &> {log.out}
        """

rule ComputeAllDeltaSVM:
    """Compute all possible and observed SVM"""
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa"
    output:
        AllSVM = pathResults + "/{TF}/deltas/ancestral_all_possible_deltaSVM.txt",
        ObsSVM = pathResults + "/{TF}/deltas/ancestral_to_observed_deltaSVM.txt"
    log: out = pathResults + "/log/{TF}/ComputeAllDeltaSVM.out"
    threads: config["nbPart"]
    params: time="2:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python  Positive_Selection_Tests/Max_LnL_Test/compute_all_deltaSVM.py {sp} \
        {sample}/{wildcards.TF} {peakType} --{cluster} -T {threads} &> {log.out}
        """

def evaluate_time(ancestral_seq):
    # Calculate time needed to perform MaxLL, considering 1500 peaks per hour and per thread.
    with open(ancestral_seq, 'r') as file:
        num_lines = sum(1 for _ in file)
        time = num_lines/ (1500*config["nbPart"])+2
        cluster_time = f"{round(time)}:00:00"
    return cluster_time

rule MaxLLTest:
    message: "Test for positive selection using Likelihood Ratio Test between 3 Maximized models"
    input:
        AllSVM = pathResults + "/{TF}/deltas/ancestral_all_possible_deltaSVM.txt",
        ObsSVM = pathResults + "/{TF}/deltas/ancestral_to_observed_deltaSVM.txt"
    output: touch(pathResults + "/{TF}/MLE_summary_{BinType}_{nbBin}bins_threshold_{threshold}.csv")
    log: out = pathResults + "/log/{TF}/MaxLLTest_{BinType}_{nbBin}_{threshold}.out"
    threads: config["nbPart"]
    params: mem="16G", threads=config["nbPart"], nbBin=config["nbBin"], BinType=config["BinType"],
            threshold=config["threshold"], time=lambda wildcards, input: evaluate_time(input.AllSVM)
    shell:
        """
        python Positive_Selection_Tests/Max_LnL_Test/MaxLL_estimation.py {sp} {sample}/{wildcards.TF} \
        --peakType {peakType} --Bins {params.BinType} --NbBin {params.nbBin} --threshold {params.threshold} \
        -T {threads} --{cluster} &> {log.out}
        """

#rule simulate_sequence:
#    """Simulate sequence evolution according to a given number of substitutions and a new optimum value"""
#    input: original_seq = "test.txt"
#    output: simulated_seq = "test.txt"
#    shell:
#        """
#        python /simulate_sequence_evolution.py {input.original_seq} {output.simulated_seq} {cluster}
#        """
