#!/usr/bin/env Rscript
library(gkmSVM)
library(BSgenome)
library(BiocManager)
args = commandArgs(trailingOnly=TRUE)

species = args[1] #"human"
sample = args[2]  #"CEBPA"
BED_file = args[3]
pathResults = args[4]

#species = "human"
#sample = "CEBPA"
#BED_file = "/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/human/human_ChIP-Seq/hsap_CEBPA.bed2"
#pathResults = "/Users/alaverre/Documents/Detecting_positive_selection/results/tools_tests/human/CEBPA/Test-models/run_15/"

if (species == "mouse"){
  genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm9.masked")
}else{
  genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38.masked")}

# For new genomes
#available.genomes()
#install("BSgenome.Mmusculus.UCSC.mm10.masked")

###############################################################################
# Generate negative sequences
genNullSeqs(BED_file, genome = genome,
            outputBedFN = paste0(pathResults,'/negSet.bed'), 
            outputPosFastaFN = paste0(pathResults,'/posSet.fa'),
            outputNegFastaFN = paste0(pathResults,'/negSet.fa'))

###############################################################################
