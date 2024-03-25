sp = config["sp"]
sample = config["sample"]
nbThreads = int(config["nbPart"])
cluster = config["cluster"]

pathResults = "./results/positive_selection/" + sp + "/" + sample
pathScripts = "./scripts/detect_positive_selection"
pathPeaks = "./results/peaks_calling/" + sp + "/" + sample

rule GenerateNegativeSeq:
    message: "Generate random sequences respecting the focal sequences composition for gkm training"
    input:
        BED_file = pathPeaks + "/{TF}.peaks.bed"
    output:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa",
        BED_UCSC = "./results/peaks_calling/" + sp + "/" + sample + "/{TF}.peaks.bed_UCSC_names"
    log: out = pathResults + "/log/{TF}/GenerateNegativeSeq.out"
    params: time="1:00:00",mem="5G",threads=1
    shell:
        """
        mkdir -p {pathResults}/{wildcards.TF}/Model/
        Rscript {pathScripts}/GenerateNegativeSeq.R {sp} {sample} {wildcards.TF} {input.BED_file} {cluster} &> {log.out}
        """

rule ModelTraining:
    message: "Training of the gkm-SVM, kmer=10, 16 threads"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.model.txt")
    log: out = pathResults + "/log/{TF}/ModelTraining.out"
    priority: 2
    params: time="15:00:00", mem="5G", threads=nbThreads
    shell:
        """
        gkmtrain -r 12 -l 10 -T {nbThreads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} &> {log.out}
        """

rule ModelValidation:
    message: "Cross-validation of the model"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.cvpred.txt")
    log: out = pathResults + "/log/{TF}/ModelValidation.out"
    params: time="48:00:00", mem="5G", threads=nbThreads
    shell:
        """
        gkmtrain -r 12 -l 10 -x 5 -T {nbThreads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} &> {log.out}
        """

rule ModelPrediction:
    message: "Generate SVM weights for all possible 10-mers"
    input:
        Model = pathResults + "/{TF}/Model/{TF}.model.txt",
        kmer_fasta = "./results/positive_selection/kmer.fa"
    output:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        Prediction_done = pathResults + "/log/{TF}/ModelPrediction_done"
    log: out = pathResults + "/log/{TF}/ModelPrediction.out"
    priority: 2
    params: time="6:00:00", mem="2G", threads=nbThreads
    shell:
        """
        gkmpredict -T {nbThreads} {input.kmer_fasta} {input.Model} {output.PredictedWeight} &> {log.out}
        """
