#!/usr/bin/env Rscript
if (!requireNamespace("gkmSVM", quietly = TRUE)) {BiocManager::install("gkmSVM", update = FALSE, ask = FALSE)}
if (!requireNamespace("BSgenome", quietly = TRUE)) {BiocManager::install("BSgenome", update = FALSE, ask = FALSE)}

suppressMessages(library(gkmSVM))
suppressMessages(library(BSgenome))
suppressMessages(library(BiocManager))
args = commandArgs(trailingOnly=TRUE)

species = args[1]
BED_file = args[2]
baseDir = args[3]
pathResults = args[4]

BED = read.table(BED_file)

# Check if chrName conversion between Ensembl and UCSC is needed
# path = ifelse(cluster=="cluster", "/work/FAC/FBM/DEE/mrobinso/evolseq/DetectPosSel", "/Users/alaverre/Documents/Detecting_positive_selection")
# pathScripts = paste0(path, "/scripts/utils/compare_genome_assemblies")

#if (species != "rat"){      #rat BED already UCSC
  #message("Convert chromosome names to UCSC convention...")
  #system(paste0("python ", pathScripts, "/convert.BED.chrNames.py ", species, " ", BED_file, " ", cluster))
#  BED_file = paste0(BED_file, "_UCSC_names")
#}

# Define mapping of species to BSgenome package names
bsgenomes <- list(
    mouse = "BSgenome.Mmusculus.UCSC.mm10.masked",
    human = "BSgenome.Hsapiens.UCSC.hg38.masked",
    dog = "BSgenome.Cfamiliaris.UCSC.canFam3.masked",
    rat = "BSgenome.Rnorvegicus.UCSC.rn6.masked",
    macaca = "BSgenome.Mmulatta.UCSC.sup2kb.rheMac8.masked",
    cat = "BSgenome.Fcatus.UCSC.sup2kb.felCat8.masked",
    cattle = "BSgenome.Btaurus.GenBank.Btau5.0.1.masked",
    rabbit = "BSgenome.Ocuniculus.UCSC.oryCun2.masked",
    chicken = "BSgenome.Ggallus.UCSC.galGal6.masked",
    spretus = "BSgenome.Mspretus.GenBank.SPRET.masked",
    caroli = "BSgenome.Mcaroli.GenBank.CAROLI.EIJ.masked",
    zebrafish = "BSgenome.Drerio.UCSC.danRer11.masked",
    drosophila = "BSgenome.Dmelanogaster.UCSC.dm6.masked"
)

pkg <- bsgenomes[[species]]

# Check if BSgenome already installed
if (!requireNamespace(pkg, quietly = TRUE)) {
    local_pkg <- file.path(baseDir, "data/BSgenome/", pkg)
    message(pkg, local_pkg)

    # Check if available in standard genomes
    if (pkg %in% available.genomes()) {BiocManager::install(pkg, update = FALSE, ask = FALSE)}

    # If local file exists, install from local path
    else if (file.exists(local_pkg)) {BiocManager::install(local_pkg, update = FALSE, ask = FALSE, type = "source")}

    else {stop("Species not recognized or BSgenome package not available locally or publicly.")}
}
genome <- getBSgenome(pkg)

###############################################################################
# Generate negative sequences
set.seed(12)

# Select at least 500 sequences or same number as positive set
max_seq = max(nrow(BED), 500)
xfold= as.integer(max_seq/nrow(BED))

genNullSeqs(BED_file, genome = genome, nMaxTrials=10, xfold=10,
            outputBedFN = paste0(pathResults,'/negSet.bed'), 
            outputPosFastaFN = paste0(pathResults,'/posSet.fa'),
            outputNegFastaFN = paste0(pathResults,'/negSet.fa'))

###############################################################################

