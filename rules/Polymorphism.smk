# Implement rules to retrieve ChiP-seq consensus peaks summits and to find their homologous using HALPER
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
    log: out = pathResults + "/log/DownloadVCF.out"
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

rule RetrieveSNPDeltaSVM:
    message: "Filter SNPs and retrieve corresponding deltaSVM and MLE estimations"
    input:
        vcf = rules.VCF_BED_overlap.output,
        AllSVM = pathResults + "/{TF}/deltas/ancestral_all_possible_deltaSVM.txt",
        ancestral_seq = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        MaxLL_estimations = pathResults + "/{TF}/MLE_summary_50bins.csv"
    output: pathPolymorphism + "/{TF}/SNP_to_deltaSVM/{chrom}.txt"
    params: time="1:00:00",mem="5G",threads=1
    shell:
        """ python peaks_evolution/SNP_to_deltaSVM.py {input.vcf} {input.AllSVM} {input.ancestral_seq} {input.MaxLL_estimations} {output} """

rule ComputeSelectionCoefficient:
    message: "Compute selection coefficient for each SNP"
    input: rules.RetrieveSNPDeltaSVM.output
    output: pathPolymorphism + "/{TF}/SelectionCoefficient/{chrom}.txt"
    params: time="1:00:00",mem="1G",threads=1
    shell:
        """ python peaks_evolution/Compute.py {input} {output} """

