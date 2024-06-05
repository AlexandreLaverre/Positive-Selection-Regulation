# Implement rules to download human polymorphism data, retrieve the estimated deltaSVM for SNP and compute their coefficient of selection
from snakemake.io import expand

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]
chroms = [f"chr{i}" for i in range(1, 22)] + ["chrX"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"
pathPolymorphism = f"../results/polymorphism_analyses/{peakType}/{sp}/{sample}"

rule DownloadVCF:
    message: "Download polymorphism data from VCF files of 1000 Genomes Project"
    output: expand("../data/polymorphism/human_1000genomes/ALL.{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", chrom=chroms)
    params: pathVCF = "../data/polymorphism/human_1000genomes"
    log: out = pathPolymorphism + "/log/DownloadVCF.out"
    shell:
        """
        mkdir -p {params.pathVCF}
        wget -A *shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz -np -nd -r -e -P {params.pathVCF} \
        robots=off http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/  &> {log.out}
        """

rule VCF_BED_overlap:
    message: "Overlap VCF with ChIP-seq peaks to filter SNP"
    input:
        vcf = '../data/polymorphism/human_1000genomes/ALL.{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz',
        BED_peaks = pathPeaks + "/save_peaks_UCSC/{TF}.peaks_UCSC_names.bed"
    output: overlap_vcf = pathPolymorphism + "/{TF}/VCF/filtered_{chrom}.vcf.gz"
    params: time="1:00:00",mem="1G",threads=1
    shell:
        """ 
        bedtools intersect -a {input.vcf} -b {input.BED_peaks} -wb -header | gzip > {output.overlap_vcf} 
        """

rule RetrieveSNPDeltaSVM_Selection:
    message: "Filter SNPs and retrieve corresponding deltaSVM and MLE estimations"
    input:
        vcf = rules.VCF_BED_overlap.output,
        AllSVM = pathResults + "/{TF}/deltas/ancestral_all_possible_deltaSVM.txt",
        focal_seq = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa",
        genome = f"../data/genome_sequences/{sp}/" + config[sp]["UCSC_Assembly"],
        MaxLL_estimations = pathResults + "/{TF}/MLE_summary_50bins.csv"
    output: pathPolymorphism + "/{TF}/SNP_to_deltaSVM/{chrom}.txt"
    log: out = pathPolymorphism + "/log/{TF}_SNP_to_delta_{chrom}.out"
    params: time="1:00:00",mem="8G",threads=1
    shell:
        """ 
        python peaks_evolution/SNP_to_deltaSVM.py {input.vcf} {input.AllSVM} {input.focal_seq} {input.genome} {input.MaxLL_estimations} {output} &> {log.out}
        """

rule MergeAllChromosome:
    message: "Compute selection coefficient for each SNP"
    input: expand(pathPolymorphism + "/{TF}/SNP_to_deltaSVM/{chrom}.txt", chrom=chroms)
    output: pathPolymorphism + "/{TF}/SNP_SelectionCoefficient.txt"
    shell:
        """cat {input} | sort -u | {output} """

