# Implement rules to simulate sequences evolution and test for positive selection:
# a comparison between the former Permutation Test and the new MaxLikelihood Ratio Test

include: '../Snakemake_rules/SVM_model.sk'

def get_paths(wildcards):
    """Input function to return all paths from the value of the wildcard 'sample'."""
    return config['data']['samples'][wildcards.sample]


rule simulate_sequence:
    """Simulate sequence evolution according to a given number of substitutions and a new optimum value"""
    input: original_seq =
    output: simulated_seq =
    shell:
        """
        python {path}/scripts/simulate_sequence_evolution.py {input.original_seq} {output.simulated_seq} {cluster}
        """

rule permutation_test:
    message: "Test for positive selection between ancestral and focal sequences using random permutations"
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        simulated_seq = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        original_seq = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa"
    output: touch(pathResults + "/{TF}/PosSelTest_deltaSVM_1000permutations.txt")
    threads: MaxThread
    shell:
        """
        python {pathScripts}/testPosSelec.py {sp} {sample} {TF} {nbRand} {cluster} --NbThread {threads}
        """

rule get_all_svm:
    """Compute all possible and observed SVM"""
    input:
        original_seq = ,
        simulated_seq =
    output:
        all_svm = ,
        obs_svm =
    threads: MaxThread
    shell:
        """
        python {path}/scripts/MaxLnl_Test/compute_all_deltaSVM.py {species} {sample} {TF} {cluster} --NbThread {threads}
        """


rule maxLL_test:
    message: "Test for positive selection using Likelihood Ratio Test between 3 Maximized models"
    input:
        all_svm = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        obs_svm = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.model.txt")
    log: out = pathResults + "/log/{TF}/ModelTraining.out"
    priority: 2
    params: time="15:00:00", mem="5G", threads=nbThreads
    shell:
        """
        gkmtrain -r 12 -l 10 -T {nbThreads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} &> {log.out}
        """
