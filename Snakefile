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
from config.config_setup import get_TFs, get_localrules, get_ruleorder
import os
configfile: 'config/TestPos.yaml'
include: 'rules/GetPeaksAlignment.smk'
include: 'rules/SVM_model.smk'
include: 'rules/FindHomologs.smk'
include: 'rules/PerformTests.smk'
include: 'rules/Polymorphism.smk'
test = get_localrules(config)
print(test)
localrules:  *get_localrules(config)
ruleorder: *get_ruleorder(config)

# Extract common variables from config
sample = config["sample"]
sp = config["sp"]
peakType = config["peakType"]
BinType = str(config["BinType"])
nbRand = str(config["nbRand"])
AncNode = str(config["AncNode"])

# Define base paths
base_results = os.path.join("..", "results")
pathResults = os.path.join(base_results, "positive_selection", peakType, sp, sample)
pathPeaks = os.path.join(base_results, "peaks_calling", peakType, sp, sample)
pathPolymorphism = os.path.join(base_results, "polymorphism_analyses", peakType, sp, sample)

# Retrieve TFs
TFs = get_TFs(config, pathPeaks)

print("Running with :", ', '.join(TFs), "transcription factors" )
########################################################################################################################
rule all:
    input :
        MaxLLTest = expand(os.path.join(pathResults, "{TF}", "Tests", f"MLE_summary_{BinType}_{AncNode}.csv"), TF=TFs),
        PosSelTest = expand(os.path.join(pathResults, "{TF}", "Tests", f"PosSelTest_deltaSVM_{nbRand}permutations_two_tailed_{AncNode}.txt") ,TF=TFs),
        archive= expand(os.path.join(pathResults, "{TF}", "alignments.archive.tar.gz"),TF=TFs)
        #model_validation = expand(pathResults + "/{TF}/Model/{TF}.cvpred.txt", TF=config["TFs"][sample]),
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
#    input:
#        PosSelTest = generate_file_paths(pathResults + "{sp}{sample}{TF}/PosSelTest_deltaSVM_" + str(config["nbRand"]) + "permutations.txt"),
#        MaxLLTest = generate_file_paths(pathResults + "{sp}{sample}{TF}/MLE_summary_" + BinType + "_" + str(nbBin) + "bins_threshold_" + str(threshold) + ".csv"),
#        archive = generate_file_paths(pathResults + "{sp}{sample}{TF}/alignments.archive.tar.gz")
