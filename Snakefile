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
from snakemake.io import expand, touch

configfile: 'config/TestEvol.yaml'
include: 'rules/GetPeaksAlignment.smk'
include: 'rules/SVM_model.smk'
include: 'rules/FindHomologs.smk'
include: 'rules/PerformTests.smk'

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"

print("Running with :", ', '.join(config["TFs"][sample]), "transcription factors" )

if cluster == "cluster":
    localrules: all, GetPeaks, BED_split, ConcatSeq, ConsensusSummits, ChromosomeCorrespondence, ConvertCoordinates
else:
    localrules: all,GetPeaks,GenerateNegativeSeq,ModelTraining,ModelValidation,ModelPrediction,BED_split,
        InferAncestralPairwise,GetSequencesMultiple,ConcatSeq,PermutationTest,ArchiveAlignments

# Define from which type of alignments ancestral sequences should be obtained
if config["AlignType"] == "pairwise":
    ruleorder: InferAncestralPairwise > GetSequencesMultiple
else:
    ruleorder: GetSequencesMultiple > InferAncestralPairwise

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
        GenomeAlignment= f"../data/genome_alignments/{sp}/triplet_ancestor.maf.gz",
        SubstiMatrixes= f"../results/substitution_matrix/{sp}/",
        ChromCorrespondence= f"../data/genome_sequences/{sp}/chromosome_correspondence.txt"
    output: Check= f"{pathResults}/log/input_check_{sp}_{sample}"
    shell:
        """
        mkdir -p {pathResults}/log
        """
