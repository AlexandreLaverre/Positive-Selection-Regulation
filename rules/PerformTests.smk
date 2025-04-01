# Implement rules to test for Positive selection and simulate sequences evolution:
from snakemake.io import touch

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]
AncNode = config["AncNode"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"

rule ComputeAllDeltaSVM:
    """Compute all possible and observed SVM"""
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_{AncNode}_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_{AncNode}_sequences_upper.fa"
    output:
        AllSVM = pathResults + "/{TF}/deltas/{AncNode}_all_possible_deltaSVM.txt",
        ObsSVM = pathResults + "/{TF}/deltas/{AncNode}_to_observed_deltaSVM.txt"
    log: out = pathResults + "/log/{TF}/ComputeAllDeltaSVM_{AncNode}.out"
    threads: config["nbPart"]
    params: time="1:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python  Positive_Selection_Tests/Max_LnL_Test/compute_all_deltaSVM.py {sp} \
        {sample}/{wildcards.TF} {peakType} --node {AncNode} --{cluster} -T {threads} > {log.out} 2>&1 
        """

def evaluate_time(ancestral_seq):
    # Calculate time needed to perform MaxLL, considering 1500 peaks per hour and per thread.
    with open(ancestral_seq, 'r') as file:
        num_lines = sum(1 for _ in file)
        time = num_lines/ (1500*config["nbPart"])+2
        cluster_time = f"{round(time)}:00:00"
    return cluster_time
#lambda wildcards, input: evaluate_time(input.AllSVM)

rule PermutationTest:
    message: "Test for positive selection between ancestral and focal sequences"
    input:
        AllSVM = pathResults + "/{TF}/deltas/{AncNode}_all_possible_deltaSVM.txt",
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_{AncNode}_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_{AncNode}_sequences_upper.fa"
    output: touch(pathResults + "/{TF}/Tests/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations_two_tailed_{AncNode}.txt")
    threads: 6 #config["nbPart"]
    priority: 2
    log: out=pathResults + "/log/{TF}/PermutationTest_{AncNode}.out"
    params: time="2:00:00", mem="5G", threads=6 #15h #config["nbPart"]
    shell:
        """
        python Positive_Selection_Tests/Permutation_Test/permutations.py {sp} {sample} {wildcards.TF} --peakType {peakType} \
        --NbRand {config[nbRand]} --node {AncNode} --{cluster} --NbThread {threads} > {log.out} 2>&1 
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
    shell:
        """
        python Positive_Selection_Tests/Max_LnL_Test/MaxLL_estimation.py {sp} {sample}/{wildcards.TF} \
        --peakType {peakType} --binType {params.BinType} --NbBin {params.nbBin} --threshold {params.threshold} \
        --node {AncNode} -T {threads} --{cluster} > {log.out} 2>&1 
        """

#rule simulate_sequence:
#    """Simulate sequence evolution according to a given number of substitutions and a new optimum value"""
#    input: original_seq = "test.txt"
#    output: simulated_seq = "test.txt"
#    shell:
#        """
#        python /simulate_sequence_evolution.py {input.original_seq} {output.simulated_seq} {cluster}
#        """
