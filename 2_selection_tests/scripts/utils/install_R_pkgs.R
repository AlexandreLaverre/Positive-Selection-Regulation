#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager", repos="https://cloud.r-project.org")}
if (!requireNamespace("gkmSVM", quietly = TRUE)) {BiocManager::install("gkmSVM", update = FALSE, ask = FALSE)}
if (!requireNamespace("BSgenome", quietly = TRUE)) {BiocManager::install("BSgenome", update = FALSE, ask = FALSE)}