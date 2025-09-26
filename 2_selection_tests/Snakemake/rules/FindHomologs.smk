# Implement rules to retrieve ChiP-seq consensus peaks summits and to find their homologous using HALPER
from snakemake.io import directory, expand
import os

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
baseDir = os.path.abspath(config["baseDir"])

pathResults = f"{baseDir}/results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"{baseDir}/results/peaks_calling/{peakType}/{sp}/{sample}"

rule ConsensusSummits:
    message: "Get consensus summits"
    input: peaks = pathPeaks + "/{TF}.peaks.bed"
    output: summits = pathPeaks + "/{TF}.consensus_summits.bed"
    log: out = pathResults + "/log/ConsensusSummits_{TF}.out"
    shell:
        """
        ../scripts/get.consensus.summits.sh {sp} {sample} > {log.out} 2>&1 
        """

rule ChromosomeCorrespondence:
    message: "Get correspondence for chromosomes names between different assemblies"
    input:
        Assembly1 = f"{baseDir}/data/genome_sequences/{sp}/" + config[sp]["Ensembl_Assembly"],
        Assembly2 = f"{baseDir}/data/genome_sequences/{sp}/" + config[sp]["UCSC_Assembly"]
    output: correspondence = f"{baseDir}/data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    log: out=pathResults + "/log/ChromosomeCorrespondence.out"
    shell:
        """
        ../scripts/utils/compare_genome_assemblies/chromosome.correspondence.sh {sp} {input.Assembly1} {input.Assembly2} Ensembl2UCSC {baseDir} > {log.out} 2>&1 
        """

rule ConvertCoordinates:
    message: "Convert coordinates to UCSC"
    input:
        peaks = pathPeaks + "/{TF}.peaks.bed",
        summits = pathPeaks + "/{TF}.consensus_summits.bed",
        correspondence = f"{baseDir}/data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    output:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    params: suffix = config[sp]["suffix"]
    shell:
        """
        python ../scripts/utils/convert.BED.chrNames.py {sp} {sample} {wildcards.TF} {params.suffix}
        """

rule runHALPER:
    message: "Get homologous peaks using HALPER"
    input:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    output: peaks = directory(baseDir +"/results/homologous_peaks/" + sp + "/{TF}/liftover/")
    log: out = pathResults + "/log/runHALPER_{TF}.out"
    container: "quay.io/comparative-genomics-toolkit/cactus:v3.0.0"
    shell:
        """
        ../scripts/peaks_evolution/run.HALPER.sh {sp} {wildcards.TF}  > {log.out} 2>&1 
        """