# Implement rules to download human polymorphism data, retrieve the estimated deltaSVM for SNP and compute their coefficient of selection
from snakemake.io import expand
import os

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
AncNode = config["AncNode"]
baseDir = os.path.abspath(config["baseDir"])

if sp == "human":
    chroms = [f"chr{i}" for i in range(1, 22)] + ["chrX"]
    vcf_prefix = "human_1000genomes/ALL."
    vcf_suffix = ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
elif sp == "drosophila":
    chroms = ["chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4"]
    vcf_prefix = "drosophila_DGRP/"
    vcf_suffix = ".dgrp2.vcf.gz"


pathResults = f"{baseDir}/results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"{baseDir}/results/peaks_calling/{peakType}/{sp}/{sample}"
pathPolymorphism = f"{baseDir}/results/polymorphism_analyses/{peakType}/{sp}/{sample}"

rule DownloadVCF:
    message: "Download polymorphism data from VCF files of 1000 Genomes Project"
    output: expand(baseDir + "/data/polymorphism/human_1000genomes/ALL.{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", chrom=chroms)
    params: pathVCF = f"{baseDir}/data/polymorphism/human_1000genomes"
    log: out = pathPolymorphism + "/log/DownloadVCF.out"
    shell:
        """
        mkdir -p {params.pathVCF}
        wget -A *shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz -np -nd -r -e -P {params.pathVCF} \
        robots=off http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/  > {log.out} 2>&1 
        """

rule VCF_BED_overlap:
    message: "Overlap VCF with ChIP-seq peaks to filter SNP"
    input:
        vcf = baseDir + '/data/polymorphism/'+ vcf_prefix + '{chrom}' + vcf_suffix,
        BED_peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed"
    output:
        overlap_vcf = pathPolymorphism + "/{TF}/VCF/filtered_{chrom}.vcf",
        overlap_vcf_gz= pathPolymorphism + "/{TF}/VCF/filtered_{chrom}.vcf.gz"
    params: time="1:00:00",mem="1G",threads=1
    log: out = pathPolymorphism + "/log/VCF_BED_{TF}_{chrom}_overlap.out"
    shell:
        """
        echo "vcf: {input.vcf}" > {log.out}
        echo "BED: {input.BED_peaks}" >> {log.out}
        mkdir -p {pathPolymorphism}/{wildcards.TF}/VCF
        bedtools intersect -a {input.vcf} -b {input.BED_peaks} -wb -header > {output.overlap_vcf} 2>> {log.out}
        gzip -c {output.overlap_vcf} > {output.overlap_vcf_gz} 2>> {log.out}
        """

rule SimpleOverlapFile:
    message: "Get unique simple file of overlapping SNP"
    input: lambda wildcards: expand(pathPolymorphism + "/{{TF}}/VCF/filtered_{chrom}.vcf.gz", chrom=chroms)
    output: pathPolymorphism + "/{TF}/VCF/overlap_peaks.txt"
    shell:
        """
        zcat {input} | grep -v '^#' | awk '{{print $NF, $2, $3}}' | sort -u > {output} 
        """

rule ComputeDeltaSVM_Reference:
    """Compute all possible SVM from reference sequence"""
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        reference_sequences = pathResults + "/{TF}/Model/posSet.fa"
    output: AllSVM = pathResults + "/{TF}/deltas/focal_ancestral_all_possible_deltaSVM.txt"
    log: out = pathResults + "/log/{TF}/ComputeDeltaSVM_Reference.out"
    threads: 2
    params: time="1:00:00", mem="5G", threads=2
    shell:
        """
        python ../scripts/RegEvol/compute_all_deltaSVM.py {sp} \
        {sample}/{wildcards.TF} {peakType} --node focal_ancestral -T {threads} > {log.out} 2>&1
        """


rule RetrieveSNPDeltaSVM_Selection:
    message: "Filter SNPs and retrieve corresponding deltaSVM and MLE estimations"
    input:
        vcf = pathPolymorphism + "/{TF}/VCF/filtered_{chrom}.vcf.gz",
        AllSVM = pathResults + "/{TF}/deltas/focal_ancestral_all_possible_deltaSVM.txt",
        focal_seq = pathResults + "/{TF}/Model/posSet.fa",
        genome = f"{baseDir}/data/genome_sequences/{sp}/" + config[sp]["UCSC_Assembly"],
        MLE = pathResults + "/{TF}/Tests/MLE_summary_exact_ranked_ancestral.csv"
    output: pathPolymorphism + "/{TF}/SNP_to_deltaSVM/{chrom}.txt"
    log: out = pathPolymorphism + "/log/{TF}_SNP_to_delta_{chrom}.out"
    params: time="1:00:00",mem="5G",threads=1
    shell:
        """ 
        python ../scripts/peaks_evolution/SNP_to_deltaSVM.py {input.vcf} {input.AllSVM} {input.focal_seq} \
        {input.genome} {input.MLE} {output} > {log.out} 2>&1 
        """

rule MergeAllChromosome:
    message: "Compute selection coefficient for each SNP"
    input: lambda wildcards: expand(pathPolymorphism + "/{{TF}}/SNP_to_deltaSVM/{chrom}.txt", chrom=chroms)
    output: pathPolymorphism + "/{TF}/SNP_SelectionCoefficient.txt"
    shell:
        """
        cat {input} | sort -u > {output}.tmp
        # ensure correct header
        (grep '^ID' {output}.tmp | head -n 1; grep -v '^ID' {output}.tmp | sort -u) > {output}
        rm {output}.tmp
        """

rule PlotSFS:
    message: "Site Frequency Spectrum"
    input:
        MLE = pathResults + "/{TF}/Tests/MLE_summary_exact_ranked_ancestral.csv",
        SelCoeff = pathPolymorphism + "/{TF}/SNP_SelectionCoefficient.txt"
    output: pathPolymorphism + "/{TF}/SFS.pdf"
    params: time="1:00:00",mem="3G",threads=1, subsample=12
    shell:
        """
        python ../scripts/peaks_evolution/plot_sfs.py --input_MLE {input.MLE} --input_SNP {input.SelCoeff} \
        --subsample {params.subsample} --output_pdf {output} 
        """