########################################################################################################################
##### PIPELINE TO DETECT POSITIVE SELECTION  #####
# Data requirements:
# input file: list of interested coordinates
# BSgenome of focal species
# pairwise or multiple genome alignments
# substitution matrix per chromosome
# chromosome correspondence if genome annotation and alignment are not from the same source

# Requirements :
# Software: PECAN + Snakemake + Exonerate + TreeTime + seqkit + lsgkm + ucsc-toupper + ucsc-MafsInRegions
# Python modules: Bio, Bio.Seq, numpy, pandas, multiprocessing.pool, alive_progress

########################################################################################################################
#! /usr/bin/env python
import glob
import os
from snakemake.io import expand, touch

configfile: 'config/TestEvol.yaml'
include: 'rules/SVM_model.smk'
inclue: 'rules/FindHomologs.smk'

sp = config["sp"]
sample = config["sample"]
AncMethod = config["AncMethod"]
cluster = config["cluster"]

suffix = "_UCSC_names" if sp in ["human", "mouse", "spretus", "caroli"] else ""
pathResults = "../results/positive_selection/NarrowPeaks/" + sp + "/" + sample
pathScripts = "../scripts/detect_positive_selection"
pathPeaks = "../results/peaks_calling/NarrowPeaks/" + sp + "/" + sample

#TFs =  list(set([os.path.basename(BED).split('_')[0] for BED in glob.glob(pathPeaks + '/bowtie2/mergedLibrary/macs2/narrowPeak/*.narrowPeak')])) ## remember to change 1 for 0
print("Running with :", ', '.join(config["TFs"][sample]), "transcription factors" )

if cluster == "cluster":
    localrules: all, GetPeaks, BED_split, ConcatSeq
else:
    localrules: all,GetPeaks,GenerateNegativeSeq,ModelTraining,ModelValidation,ModelPrediction,BED_split,
        InferAncestralPairwise,GetSequencesMultiple,ConcatSeq,TestPosSel,ArchiveAlignments

# Define from which type of alignments ancestral sequences should be obtained
if config["AlignType"] == "pairwise":
    ruleorder: InferAncestralPairwise > GetSequencesMultiple
else:
    ruleorder: GetSequencesMultiple > InferAncestralPairwise

workdir: pathResults
########################################################################################################################

rule all:
    input :
        PosSelTest = expand(pathResults + "/{TF}/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations.txt", TF=config["TFs"][sample]),
        archive= expand(pathResults + "/{TF}/alignments.archive.tar.gz", TF=config["TFs"][sample]),
        #model_validation = expand(pathResults + "/{TF}/Model/{TF}.cvpred.txt", TF=config["TFs"][sample])

rule check_input_data:
    message: "Check if all the required data are present before starting."
    input:
        PeaksFolder = f"{pathPeaks}/bowtie2/mergedLibrary/macs2/narrowPeak/",
        GenomeAlignment=f"data/genome_alignments/{sp}/triplet_ancestor.maf.gz",
        SubstiMatrixes=f"results/substitution_matrix/{sp}/",
        ChromCorrespondence=f"data/genome_sequences/{sp}/chromosome_correspondence.txt"
    output: Check=f"{pathResults}/log/input_check_{sp}_{sample}"
    shell:
        """
        mkdir -p {pathResults}/log
        """

rule GetPeaks:
    message: "Retrieve ChIP peaks with a meaningful ID"
    input: PeaksFolder = f"{pathPeaks}/bowtie2/mergedLibrary/macs2/narrowPeak/",
           InputCheck = f"{pathResults}/log/input_check_{sp}_{sample}"
    output: Peaks = pathPeaks + "/{TF}.peaks.bed"
    shell:
        """
        if [ -d {input.PeaksFolder}/consensus/{wildcards.TF} ]; then
            # peaks come from a consensus of several samples
	        cp {input.PeaksFolder}/consensus/{wildcards.TF}/{wildcards.TF}.consensus_peaks.bed {output.Peaks} 
        else
            # peaks come from an unique sample
            cp {input.PeaksFolder}/{wildcards.TF}*.narrowPeak {output.Peaks}
        fi
        
        # Add meaningful ID for each peak
        cut -f 1-4 {output.Peaks} | sed 's/\t/:/g' > {pathPeaks}/{wildcards.TF}_IDs
        cut -f 1-3 {output.Peaks} > {pathPeaks}/{wildcards.TF}_coord
        paste {pathPeaks}/{wildcards.TF}_coord {pathPeaks}/{wildcards.TF}_IDs > {output.Peaks}
        rm {pathPeaks}/{wildcards.TF}_coord {pathPeaks}/{wildcards.TF}_IDs
        """

rule ConvertCoordinates:
    message: "Convert coordinates to UCSC for human and mice"
    input:
        peaks = pathPeaks + "/{TF}.peaks.bed",
        summits = pathPeaks + "/consensus/{TF}/{TF}.consensus_summits.bed",
        correspondence = f"data/genome_sequences/{sp}/chromosome_correspondence_Ensembl2UCSC.txt"
    output:
        peaks = pathPeaks + "/{TF}.peaks_UCSC_names.bed",
        summits = pathPeaks + "/{TF}.consensus_summits_UCSC_names.bed"
    shell:
        """
        python scripts/utils/convert.BED.chrNames.py {sp} {sample} {cluster}
        """

rule BED_split:
    message: "Split the list of coordinates for parallelization"
    input: BED_file = pathPeaks + "/{TF}.peaks" + suffix + ".bed"
    output: touch(pathResults + "/log/{TF}/part{part}")
    shell:
        """
        mkdir -p {pathResults}/log/
        {config[split][cluster]} -d -n l/{config[nbPart]} {input.BED_file} {pathResults}/log/{wildcards.TF}/part1
        """

rule InferAncestralPairwise:
    message: "Infer ancestral sequences from pairwise alignments"
    input:
        BED_file_part = pathResults + "/log/{TF}/part{part}",
        Positive_seq = pathResults + "/{TF}/Model/posSet.fa"
    output: touch(pathResults + "/log/{TF}/GetAncestral_part{part}_done")
    log: out = pathResults + "/log/{TF}/GetAncestral_part{part}.out"
    params: time="2:00:00",mem="1G",threads=1
    shell:
        """
        mkdir -p {pathResults}/{wildcards.TF}/Alignments/
        python {pathScripts}/InferAncestralPairwise.py {sp} {sample} {wildcards.TF} {input.BED_file_part} {AncMethod} {cluster} &> {log.out}
        """

rule GetSequencesMultiple:
    message: "Retrieve focal and ancestral sequences from multiple whole-genome alignment"
    input:
        BED_file_part = pathResults + "/log/{TF}/part{part}"
    output: touch(pathResults + "/log/{TF}/GetAncestral_part{part}_done")
    log: out = pathResults + "/log/{TF}/extract_sequences_from_MAF_part{part}.out"
    params: time="2:00:00",mem="1G",threads=1
    shell:
        """
        {pathScripts}/extract_sequences_from_MAF.sh {sp} {sample} {wildcards.TF} {input.BED_file_part} {cluster} &> {log.out}
        """

rule ConcatSeq:
    message: "Concatenate and sort all coordinates"
    input:
        Ancestral= lambda wildcards: expand(pathResults + "/log/{TF}/GetAncestral_part{part}_done",part=range(100,100 + config["nbPart"]),TF=wildcards.TF)
    output:
        concat_ancestral = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        concat_focal = pathResults + "/{TF}/sequences/all_focal_sequences.fa",
        concat_focal_filtered = pathResults + "/{TF}/sequences/filtered_focal_sequences.fa",
        concat_focal_upper = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa",
        concat_sister= pathResults + "/{TF}/sequences/all_sister_sequences.fa",
        concat_sister_filtered= pathResults + "/{TF}/sequences/filtered_sister_sequences.fa",
        concat_sister_upper= pathResults + "/{TF}/sequences/filtered_sister_sequences_upper.fa",
    params: time="1:00:00",mem="1G",threads=1
    shell:
        """
        cd {pathResults}/{wildcards.TF}/Alignments/
        mkdir -p {pathResults}/{wildcards.TF}/sequences/
        # Get all non empty ancestral sequences in one file
        find ancestral_sequences -name "*nogap.fa" -size +0 | xargs basename -s _nogap.fa > list_ancestral.txt
        find ancestral_sequences -name "*nogap.fa" -size +0 -exec cat {{}} + > {output.concat_ancestral}

        # Get all corresponding focal sequences in one file
        find focal_sequences -name "*nogap.fa" -size +0 -exec cat {{}} + > {output.concat_focal}
        seqtk subseq {output.concat_focal} {pathResults}/{wildcards.TF}/Alignments/list_ancestral.txt > {output.concat_focal_filtered}
        
        # Get all corresponding sister species's sequences in one file
        find sister_sequences -name "*nogap.fa" -size +0 -exec cat {{}} + > {output.concat_sister}
        seqtk subseq {output.concat_sister} {pathResults}/{wildcards.TF}/Alignments/list_ancestral.txt > {output.concat_sister_filtered}
        
        # Make sequences in uppercase to remove potential soft repeat mask 
        awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {output.concat_focal_filtered} > {output.concat_focal_upper}
        awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' {output.concat_sister_filtered} > {output.concat_sister_upper}
        """

rule ArchiveAlignments:
    message: "Create an archive file containing all alignments"
    input: focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa"
    output: archive = pathResults + "/{TF}/alignments.archive.tar.gz"
    params: time="1:00:00",mem="1G",threads=1
    shell:
        """
        tar -czf {output.archive} -C {pathResults}/{wildcards.TF} Alignments/
        rm -r {pathResults}/{wildcards.TF}/Alignments/
        """

rule TestPosSel:
    message: "Test for positive selection between ancestral and focal sequences"
    input:
        PredictedWeight = pathResults + "/{TF}/Model/kmer_predicted_weight.txt",
        ancestral_sequences = pathResults + "/{TF}/sequences/filtered_ancestral_sequences.fa",
        focal_sequences = pathResults + "/{TF}/sequences/filtered_focal_sequences_upper.fa"
    output: touch(pathResults + "/{TF}/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations.txt")
    threads: config["nbPart"]
    priority: 2
    params: time="15:00:00", mem="5G", threads=config["nbPart"]
    shell:
        """
        python {pathScripts}/testPosSelec.py {sp} {sample} {wildcards.TF} {config[nbRand]} {cluster} --NbThread {threads}
        """

