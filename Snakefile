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
include: 'rules/Polymorphism.smk'

sp = config["sp"]
sample = config["sample"]
peakType = config["peakType"]
cluster = config["cluster"]
threshold = str(config["threshold"])

pathResults = f"../results/positive_selection/{peakType}/{sp}/{sample}"
pathPeaks = f"../results/peaks_calling/{peakType}/{sp}/{sample}"
pathPolymorphism = f"../results/polymorphism_analyses/{peakType}/{sp}/{sample}"

print("Running with :", ', '.join(config["TFs"][sample]), "transcription factors" )

if cluster == "cluster":
    localrules: all, GetPeaks, BED_split, ConcatSeq, ConsensusSummits, ModelPrediction, ChromosomeCorrespondence, ConvertCoordinates, DownloadVCF, MergeAllChromosome
else:
    localrules: all, GetPeaks,GenerateNegativeSeq,ModelTraining,ModelValidation,ModelPrediction,BED_split,
        InferAncestralPairwise,GetSequencesMultiple,ConcatSeq,PermutationTest,ArchiveAlignments, MergeAllChromosome

# Define from which type of alignments ancestral sequences should be obtained
if config["AlignType"] == "pairwise":
    ruleorder: InferAncestralPairwise > GetSequencesMultiple
else:
    ruleorder: GetSequencesMultiple > InferAncestralPairwise

########################################################################################################################
rule all:
    input :
        MaxLLTest = expand(pathResults + "/{TF}/Tests/MLE_summary_" + str(config["BinType"]) + "_" + str(config["nbBin"]) + "bins_threshold_" + str(config["threshold"]) + "_last.csv", TF=config["TFs"][sample]),
        PosSelTest = expand(pathResults + "/{TF}/Tests/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations_last.txt",TF=config["TFs"][sample]),
        archive= expand(pathResults + "/{TF}/alignments.archive.tar.gz",TF=config["TFs"][sample]),

        model_validation = expand(pathResults + "/{TF}/Model/{TF}.cvpred.txt", TF=config["TFs"][sample]),
        #SNP_to_delta= expand(pathPolymorphism + "/{TF}/SNP_SelectionCoefficient.txt", TF=config["TFs"][sample])

########################################################################################################################
# TMP run on all species
#def generate_file_paths(file_pattern):
#    return [
#        file_pattern.format(sp=sp, sample=sample, TF=TF)
#        for sp in config["species"]
#        for sample in config[sp]["sample"]
#        for TF in config["TFs"][sample]
#    ]

#rule all:
    input:
        PosSelTest = generate_file_paths(pathResults + "{sp}{sample}{TF}/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations.txt"),
        MaxLLTest = generate_file_paths(pathResults + "{sp}{sample}{TF}/MLE_summary_" + BinType + "_" + str(nbBin) + "bins_threshold_" + str(threshold) + ".csv"),
        archive = generate_file_paths(pathResults + "{sp}{sample}{TF}/alignments.archive.tar.gz")
