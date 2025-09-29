#!/usr/bin/env Rscript
if (!requireNamespace("gkmSVM", quietly = TRUE)) {BiocManager::install("gkmSVM", update = FALSE, ask = FALSE)}
if (!requireNamespace("BSgenome", quietly = TRUE)) {BiocManager::install("BSgenome", update = FALSE, ask = FALSE)}

suppressMessages(library(gkmSVM))
suppressMessages(library(BSgenome))
args = commandArgs(trailingOnly=TRUE)

species = args[1]
BED_file = args[2]
baseDir = args[3]
pathResults = args[4]

BED = read.table(BED_file)

# Define mapping of species to BSgenome package names
bsgenomes <- list(
    mouse = "BSgenome.Mmusculus.UCSC.mm10",
    human = "BSgenome.Hsapiens.UCSC.hg38",
    dog = "BSgenome.Cfamiliaris.UCSC.canFam3",
    rat = "BSgenome.Rnorvegicus.UCSC.rn6",
    macaca = "BSgenome.Mmulatta.UCSC.sup2kb.rheMac8",
    cat = "BSgenome.Fcatus.UCSC.sup2kb.felCat8",
    cattle = "BSgenome.Btaurus.GenBank.Btau5.0.1",
    rabbit = "BSgenome.Ocuniculus.UCSC.oryCun2",
    chicken = "BSgenome.Ggallus.UCSC.galGal6",
    spretus = "BSgenome.Mspretus.GenBank.SPRET",
    caroli = "BSgenome.Mcaroli.GenBank.CAROLI.EIJ",
    zebrafish = "BSgenome.Drerio.UCSC.danRer11",
    drosophila = "BSgenome.Dmelanogaster.UCSC.dm6"
)

pkg <- bsgenomes[[species]]
pkg_masked <- paste0(pkg, ".masked")

# Check if BSgenome already installed
for (p in c(pkg, pkg_masked)) {
    local_pkg <- file.path(baseDir, "data/BSgenome", paste0(p, "_1.0.tar.gz"))
    public_genomes <- available.genomes()
    if (!requireNamespace(p, quietly = TRUE)) {
        if (p %in% public_genomes) {
            message("Installing ", p, " from Bioconductor...")
            BiocManager::install(p, update = FALSE, ask = FALSE)
        } else if (file.exists(local_pkg)) {
            message("Installing ", p, " from local file: ", local_pkg)
            install.packages(local_pkg, repos = NULL, type = "source")
        } else {
            stop("Package ", p, " not found locally or in Bioconductor.")
        }
    } else {
        message(p, " already installed.")
    }
}

genome <- getBSgenome(pkg_masked)

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

