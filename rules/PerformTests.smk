# Implement rules to simulate sequences evolution and test for positive selection:
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
    params: time="15:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python Positive_Selection_Tests/Permutation_Test/permutations.py {sp} {sample} {wildcards.TF} {config[nbRand]} {cluster} --NbThread {threads}
        """

#rule simulate_sequence:
#    """Simulate sequence evolution according to a given number of substitutions and a new optimum value"""
#    input: original_seq = "test.txt"
#    output: simulated_seq = "test.txt"
#    shell:
#        """
#        python /simulate_sequence_evolution.py {input.original_seq} {output.simulated_seq} {cluster}
#        """


#rule get_all_svm:
#    """Compute all possible and observed SVM"""
#    input:
#        original_seq = "test.txt",
#        simulated_seq ="test.txt"
#    output:
#        all_svm = "test.txt",
#        obs_svm = "test.txt"
#    threads: config["nbPart"]
#    shell:
#        """
#        python {path}/scripts/MaxLnl_Test/compute_all_deltaSVM.py {sp} {sample} {wildcards.TF} {cluster} --NbThread {threads}
#        """


#rule maxLL_test:
#    message: "Test for positive selection using Likelihood Ratio Test between 3 Maximized models"
#    input:
#        all_svm = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
#        obs_svm = pathResults + "/{TF}/Model/negSet.fa"
#    output: touch(pathResults + "/{TF}/Model/{TF}.model.txt")
#    log: out = pathResults + "/log/{TF}/ModelTraining.out"
#    priority: 2
#    params: time="15:00:00", mem="5G", threads=config["nbPart"]
#    shell:
#        """
#        gkmtrain -r 12 -l 10 -T {nbThreads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} &> {log.out}
#        """
