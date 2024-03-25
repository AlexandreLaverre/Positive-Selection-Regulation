sp = config["sp"]
sample = config["sample"]
cluster = config["cluster"]

pathResults = "../results/positive_selection/NarrowPeaks/" + sp + "/" + sample
pathPeaks = "../results/peaks_calling/NarrowPeaks/" + sp + "/" + sample

rule ConsensusSummits:
    message: "Get consensus summits"
    input: peaks = pathPeaks + "/consensus/{TF}/{TF}.consensus_peaks.bed"
    output: summits = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed"
    log: out = pathResults + "/log/ConsensusSummits_{TF}.out"
    shell:
        """
        /ChIPseq_analyses/get.consensus.summits.sh {sp} {sample} {cluster} &> {log.out}
        """

rule ChromosomeCorrespondence:
    message: "Get correspondence for chromosomes names between different assemblies"
    input:
        Assembly1 = f"../data/genome_sequences/{sp}/" + config[sp]["Ensembl_Assembly"],
        Assembly2 = f"../data/genome_sequences/{sp}/" + config[sp]["UCSC_Assembly"]
    output: correspondence = f"../data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    shell:
        """
        /utils/compare_genome_assemblies/chromosome.correspondence.sh {sp} {input.Assembly1} {input.Assembly2} Ensembl2UCSC {cluster}
        """

rule ConvertCoordinates:
    message: "Convert coordinates to UCSC for human and mice"
    input:
        peaks = pathPeaks + "/{TF}.peaks.bed",
        summits = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed",
        correspondence = f"../data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    output:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    shell:
        """
        python scripts/utils/convert.BED.chrNames.py {sp} {sample} {cluster}
        """

rule runHALPER:
    message: "Get consensus summits"
    input: BED_file = pathPeaks + "/{TF}.peaks.bed"
    output: peaks = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed"
    log: out = pathResults + "/log/runHALPER_{TF}.out"
    shell:
        """
        /peaks_evolution/run.HALPER.sh {sp} {wildcards.TF} {cluster} &> {log.out}
        """