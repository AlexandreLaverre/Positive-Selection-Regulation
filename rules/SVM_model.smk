from snakemake.io import touch, directory

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
nbThreads = config["nbPart"]
cluster = config["cluster"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"
pathScripts = "../scripts/Positive_Selection_Tests/"


rule GenerateNegativeSeq:
    message: "Generate random sequences respecting the focal sequences composition for gkm training"
    input:
        BED = pathPeaks + "/{TF}.peaks_UCSC_names.bed"
    output:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    log: out = pathResults + "/log/{TF}/GenerateNegativeSeq.out"
    priority: 2
    params: time="1:00:00",mem="10G",threads=1 #5h
    shell:
        """
        pathModel="{pathResults}/{wildcards.TF}/Model/"
        mkdir -p $pathModel
        Rscript {pathScripts}/generate_negative_sequence.R {sp} {input.BED} $pathModel {cluster} > {log.out} 2>&1 
        """

rule ModelTraining:
    message: "Training of the gkm-SVM, kmer=10"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.model.txt")
    log: out = pathResults + "/log/{TF}/ModelTraining.out"
    priority: 2
    threads: config["ModelThreads"]
    params: time="2:00:00", mem="5G", threads=config["ModelThreads"] #24h
    shell:
        """
        gkmtrain -r 12 -l 10 -T {threads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} > {log.out} 2>&1 
        """

rule ModelValidation:
    message: "Cross-validation of the model"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.cvpred.txt")
    log: out = pathResults + "/log/{TF}/ModelValidation.out"
    threads: config["ModelThreads"]
    params: time="48:00:00", mem="5G", threads=config["ModelThreads"]
    shell:
        """
        gkmtrain -r 12 -l 10 -x 5 -T {threads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} > {log.out} 2>&1 
        """

rule ModelPrediction:
    message: "Generate SVM weights for all possible 10-mers, 4 threads"
    input:
        Model = pathResults + "/{TF}/Model/{TF}.model.txt",
        kmer_fasta = "../results/positive_selection/kmer.fa"
    output: PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt"
    log: out = pathResults + "/log/{TF}/ModelPrediction.out"
    priority: 2
    params: time="1:00:00", mem="2G", threads=1
    shell:
        """
        gkmpredict -T 1 {input.kmer_fasta} {input.Model} {output.PredictedWeight} > {log.out} 2>&1 
        """
