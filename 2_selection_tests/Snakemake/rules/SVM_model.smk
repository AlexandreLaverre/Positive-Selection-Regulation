from snakemake.io import touch, directory
import os

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
AncNode = config["AncNode"]
baseDir = os.path.abspath(config["baseDir"])

pathResults = f"{baseDir}/results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"{baseDir}/results/peaks_calling/{peakType}/{sp}/{sample}"

rule install_r_pkgs:
    output: touch(pathResults + "/log/r_pkgs_installed")
    conda: "../../envs/training_gkm.yaml"
    shell:
        """
        Rscript scripts/install_R_pkgs.R
        touch {output}
        """

rule GenerateNegativeSeq:
    message: "Generate random sequences respecting the focal sequences composition for gkm training"
    input:
        BED = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        pkg_install = pathResults + "/log/r_pkgs_installed"
    output:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    log: out = pathResults + "/log/{TF}/GenerateNegativeSeq.out"
    priority: 10
    params: time="5:00:00",mem="10G",threads=1
    conda: "../../envs/training_gkm.yaml"
    shell:
        """
        pathModel="{pathResults}/{wildcards.TF}/Model/"
        mkdir -p $pathModel
        Rscript ../scripts/generate_negative_sequence.R {sp} {input.BED} $pathModel > {log.out} 2>&1 
        """

rule ModelTraining:
    message: "Training of the gkm-SVM, kmer=10"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.model.txt")
    log: out = pathResults + "/log/{TF}/ModelTraining.out"
    priority: 10
    threads: config["ModelThreads"]
    params: time="24:00:00", mem="5G", threads=config["ModelThreads"]
    conda: "../../envs/training_gkm.yaml"
    shell:
        """
        gkmtrain -r 12 -l 10 -T {threads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} > {log.out} 2>&1 || exit 1
        """

rule ModelValidation:
    message: "Cross-validation of the model"
    input:
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa",
        Negative_seq = pathResults + "/{TF}/Model/negSet.fa"
    output: touch(pathResults + "/{TF}/Model/{TF}.cvpred.txt")
    log: out = pathResults + "/log/{TF}/ModelValidation.out"
    threads: config["ModelThreads"]
    params: time="20:00:00", mem="5G", threads=config["ModelThreads"]
    conda: "../../envs/training_gkm.yaml"
    shell:
        """
        gkmtrain -r 12 -l 10 -x 5 -T {threads} {input.Positive_seq} {input.Negative_seq} {pathResults}/{wildcards.TF}/Model/{wildcards.TF} > {log.out} 2>&1 || exit 1
        """

rule ModelPrediction:
    message: "Generate SVM weights for all possible 10-mers, 4 threads"
    input:
        Model = pathResults + "/{TF}/Model/{TF}.model.txt",
        kmer_fasta = baseDir + "/results/positive_selection/kmer.fa"
    output: PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt"
    log: out = pathResults + "/log/{TF}/ModelPrediction.out"
    priority: 2
    params: time="1:00:00", mem="2G", threads=1
    conda: "../../envs/training_gkm.yaml"
    shell:
        """
        gkmpredict -T 1 {input.kmer_fasta} {input.Model} {output.PredictedWeight} > {log.out} 2>&1 
        """
