#!/usr/bin/env Rscript
suppressMessages(library(gkmSVM))
suppressMessages(library(BSgenome))
suppressMessages(library(BiocManager))
#library(reticulate)
#use_python("/Users/alaverre/miniconda3/bin/python")
args = commandArgs(trailingOnly=TRUE)

species = args[1] #"human"
sample = args[2]  #"CEBPA"
BED_file = args[3]
pathResults = args[4]
pathScripts = "/Users/alaverre/Documents/Detecting_positive_selection/scripts/compare_genome_assemblies"

# Check if chrName conversion between Ensembl and UCSC is needed
BED = read.table(BED_file)
if ("Ensembl" %in% seqlevelsStyle(BED$V1)){
  # Run chromosome names correspondence
  system(paste0("python ", pathScripts, "/convert.BED.chrNames.py ", species, " ", BED_file))
  BED_file = paste0(BED_file, "_UCSC_names")
}

if (species == "mouse"){genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm9.masked")}
if (species == "human"){genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38.masked")}
if (species == "dog"){genome <- getBSgenome("BSgenome.Cfamiliaris.UCSC.canFam3.masked")}

# For new genomes
#available.genomes()
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")

###############################################################################
# Generate negative sequences
genNullSeqs(BED_file, genome = genome,
            outputBedFN = paste0(pathResults,'/negSet.bed'), 
            outputPosFastaFN = paste0(pathResults,'/posSet.fa'),
            outputNegFastaFN = paste0(pathResults,'/negSet.fa'))

###############################################################################
