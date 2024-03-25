from snakemake.io import expand

sp = config["sp"]
sample = config["sample"]
cluster = config["cluster"]

pathResults = "results/positive_selection/NarrowPeaks/" + sp + "/" + sample
pathPeaks = "results/peaks_calling/NarrowPeaks/" + sp + "/" + sample

rule ConsensusSummits:
    message: "Get consensus summits"
    input: peaks = pathPeaks + "/consensus/{TF}/{TF}.consensus_peaks.bed"
    output: summits = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed"
    log: out = pathResults + "/log/ConsensusSummits_{TF}.out"
    shell:
        """
        scripts/ChIPseq_analyses/get.consensus.summits.sh {sp} {sample} {cluster} &> {log.out}
        """

rule ChromosomeCorrespondence:
    message: "Get correspondence for chromosomes names between different assemblies"
    input:
        Assembly1 = config[sp]["Ensembl_Assembly"],
        Assembly2 = config[sp]["UCSC_Assembly"]
    output: correspondence = f"data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    shell:
        """
        scripts/utils/compare_genome_assemblies/chromosome.correspondence.sh {sp} {input.Assembly1} {input.Assembly2} Ensembl2UCSC {cluster}
        """


rule runHALPER:
    message: "Get consensus summits"
    input: BED_file = pathPeaks + "/{TF}.peaks.bed"
    output: peaks = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed"
    log: out = pathResults + "/log/runHALPER_{TF}.out"
    shell:
        """
        scripts/peaks_evolution/run.HALPER.sh {sp} {wildcards.TF} {cluster} &> {log.out}
        """