#!/usr/bin/env Rscript
suppressMessages(library(gkmSVM))
suppressMessages(library(BSgenome))
suppressMessages(library(BiocManager))
args = commandArgs(trailingOnly=TRUE)

species = args[1] #"human"
sample = args[2]  #"Wilson"
TF = args[3]  #"CEBPA"
BED_file = args[4]
cluster = args[5]

path = ifelse(cluster=="cluster", "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel", "/Users/alaverre/Documents/Detecting_positive_selection")

pathScripts = paste0(path, "/scripts/compare_genome_assemblies")
pathResults = paste0(path, "/results/positive_selection/", species, "/", sample, "/", TF, "/Model")

BED = read.table(BED_file)

# Check if chrName conversion between Ensembl and UCSC is needed
#if (species != "rat"){      #rat BED already UCSC
#  message("Convert chromosome names to UCSC convention...")
#  system(paste0("python ", pathScripts, "/convert.BED.chrNames.py ", species, " ", BED_file, " ", cluster))
#  BED_file = paste0(BED_file, "_UCSC_names")
#}

if (species == "mouse"){genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")}
if (species == "human"){genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38.masked")}
if (species == "dog"){genome <- getBSgenome("BSgenome.Cfamiliaris.UCSC.canFam3.masked")}
if (species == "rat"){genome <- getBSgenome("BSgenome.Rnorvegicus.UCSC.rn6.masked")}
if (species == "macaca"){genome <- getBSgenome("BSgenome.Mmulatta.UCSC.sup2kb.rheMac8.masked")}
if (species == "cat"){genome <- getBSgenome("BSgenome.Fcatus.UCSC.sup2kb.felCat8.masked")}
if (species == "cattle"){genome <- getBSgenome("BSgenome.Btaurus.GenBank.Btau5.0.1.masked")}
if (species == "rabbit"){genome <- getBSgenome("BSgenome.Ocuniculus.UCSC.oryCun2.masked")}
if (species == "chicken"){genome <- getBSgenome("BSgenome.Ggallus.UCSC.galGal6.masked")}
#if (species == "pig"){genome <- getBSgenome("BSgenome.Sscrofa.UCSC.susScr3.masked")}
if (species == "spretus"){genome <- getBSgenome("BSgenome.Mspretus.GenBank.SPRET.masked")}
if (species == "caroli"){genome <- getBSgenome("BSgenome.Mcaroli.GenBank.CAROLI.EIJ.masked")}


# For new genomes, check if available first, else create it via BSgenome forge
#available.genomes()
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")

###############################################################################
# Generate negative sequences
set.seed(12)

genNullSeqs(BED_file, genome = genome, nMaxTrials=100,
            outputBedFN = paste0(pathResults,'/negSet.bed'), 
            outputPosFastaFN = paste0(pathResults,'/posSet.fa'),
            outputNegFastaFN = paste0(pathResults,'/negSet.fa'))

###############################################################################
