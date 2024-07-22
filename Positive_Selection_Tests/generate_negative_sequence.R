#!/usr/bin/env Rscript
suppressMessages(library(gkmSVM))
suppressMessages(library(BSgenome))
suppressMessages(library(BiocManager))
args = commandArgs(trailingOnly=TRUE)

species = args[1] #"human"
BED_file = args[2]
pathResults = args[3]
cluster = args[4]

BED = read.table(BED_file)

# Check if chrName conversion between Ensembl and UCSC is needed
# path = ifelse(cluster=="cluster", "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel", "/Users/alaverre/Documents/Detecting_positive_selection")
# pathScripts = paste0(path, "/scripts/utils/compare_genome_assemblies")

#if (species != "rat"){      #rat BED already UCSC
  #message("Convert chromosome names to UCSC convention...")
  #system(paste0("python ", pathScripts, "/convert.BED.chrNames.py ", species, " ", BED_file, " ", cluster))
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
if (species == "zebrafish"){genome <- getBSgenome("BSgenome.Drerio.UCSC.danRer11.masked")}
if (species == "drosophila"){genome <- getBSgenome("BSgenome.Dmelanogaster.UCSC.dm6.masked")}

#seqnames(genome)
# For new genomes, check if available first, else create it via BSgenome forge
#available.genomes()
#BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6.masked")

###############################################################################
# Generate negative sequences
set.seed(12)

# Select at least 5000 sequences or same number as positive set
max_seq = max(nrow(BED), 5000)
xfold= as.integer(max_seq/nrow(BED))

genNullSeqs(BED_file, genome = genome, nMaxTrials=50, xfold=xfold,
            outputBedFN = paste0(pathResults,'/negSet.bed'), 
            outputPosFastaFN = paste0(pathResults,'/posSet.fa'),
            outputNegFastaFN = paste0(pathResults,'/negSet.fa'))

###############################################################################
