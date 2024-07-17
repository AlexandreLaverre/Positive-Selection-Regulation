# Implement rules to retrieve ChiP-seq consensus peaks summits and to find their homologous using HALPER
from snakemake.io import directory, expand

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"


rule ConsensusSummits:
    message: "Get consensus summits"
    input: peaks = expand(pathPeaks + "/{TF}.peaks.bed", TF=config["TFs"][sample])
    output: summits = expand(pathPeaks + "/{TF}.consensus_summits.bed", TF=config["TFs"][sample])
    log: out = pathResults + "/log/ConsensusSummits.out"
    shell:
        """
        ./ChIPseq_analyses/get.consensus.summits.sh {sp} {sample} {cluster} &> {log.out}
        """

rule ChromosomeCorrespondence:
    message: "Get correspondence for chromosomes names between different assemblies"
    input:
        Assembly1 = f"../data/genome_sequences/{sp}/" + config[sp]["Ensembl_Assembly"],
        Assembly2 = f"../data/genome_sequences/{sp}/" + config[sp]["UCSC_Assembly"]
    output: correspondence = f"../data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    log: out=pathResults + "/log/ChromosomeCorrespondence.out"
    shell:
        """
        ./utils/compare_genome_assemblies/chromosome.correspondence.sh {sp} {input.Assembly1} {input.Assembly2} Ensembl2UCSC {cluster} &> {log.out}
        """

rule ConvertCoordinates:
    message: "Convert coordinates to UCSC"
    input:
        peaks = pathPeaks + "/{TF}.peaks.bed",
        summits = pathPeaks + "/{TF}.consensus_summits.bed",
        correspondence = f"../data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    output:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    shell:
        """
        python utils/convert.BED.chrNames.py {sp} {sample} {wildcards.TF} {cluster}
        """

rule runHALPER:
    message: "Get consensus summits"
    input:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    output: peaks = directory("../results/homologous_peaks/" + sp + "/{TF}/liftover/")
    log: out = pathResults + "/log/runHALPER_{TF}.out"
    shell:
        """
        peaks_evolution/run.HALPER.sh {sp} {wildcards.TF} {cluster} &> {log.out}
        """