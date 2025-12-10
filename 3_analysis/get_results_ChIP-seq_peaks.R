library(stringr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Biostrings)
library(qvalue)
library(gap)

options(scipen = 1)
path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

species=c("caroli", "spretus", "dog", "cat", "human", "macaca", "rabbit", "rat", "mouse", "chicken", "drosophila")
sample=c("Wilson", "Myers", "Rensch", "Schmidt10", "Schmidt12", "Stefflova", "Ni12")
TFs=c("CEBPA", "FOXA1", "HNF4A", "HNF6", "CTCF")
maxSub=150
minSub=2
pval_threshold=0.01
fdr_threshold=0.05

method="exact_ranked" # quantile_50bins_threshold_0.01
################################################################################
# Get an object of all merged results

AllMerged = paste0(path, "allMLE_list_", method, "_", minSub, "sub.Rds")
if (!file.exists(AllMerged)){
  all_MLE_list <- list()
  for (sp in species){
      for (samp in sample){
        for (TF in TFs){
          MLE.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/MLE_summary_", method, ".csv")
          delta.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/deltas/ancestral_to_observed_deltaSVM.txt")
          permut.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
          stats.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/sequences/focal_substitutions_stats.txt")
          
          if (!file.exists(MLE.file)){next}
          print(paste(sp, samp, TF))
          
          MLE <- read.csv(MLE.file, h=T, row.names = 1)
          
          # Stats substitutions
          seq_stats <- read.table(stats.file, h=T, row.names = 1)
          seq_stats <- seq_stats[row.names(MLE),]
          
          # Observed deltas
          obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
          deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
          row.names(deltas) <- deltas$seq_name
          deltas <- deltas[row.names(MLE),]
          
          MLE$SVM <- deltas$SVM
          MLE$deltaSVM <- deltas$deltaSVM
          names <- gsub("_.*:\\d+:\\d+:", "_", rownames(MLE))
          MLE$ID <- gsub("[_\\:]([[:alpha:]]).*_.*", "", names)
          MLE$peaks_ID <- gsub(".*:\\d+:\\d+[_\\:]", "", rownames(MLE))
          
          # Filtered out according to indel and Nsub
          MLE$original_length <- as.numeric(str_split_i(MLE$ID, ":", 3))-as.numeric(str_split_i(MLE$ID, ":", 2))
          MLE$indel <- MLE$original_length - seq_stats$Length
          
          MLE$Diffdelta <- MLE$deltaSVM-MLE$SumObs
          MLE$Prop.indel <- (MLE$indel*100)/MLE$original_length
          MLE$original_sample_size <- nrow(MLE)
          MLE <- MLE[which(MLE$Prop.indel<10),]
          MLE <- MLE[which(MLE$Nmut >= minSub),]
          
          # FDR corrections
          MLE$FDR_null_pos <- p.adjust(MLE$p_value_null_pos, method="fdr")
          MLE$FDR_null_purif <- p.adjust(MLE$p_value_null_purif, method="fdr")
          MLE$FDR_purif_pos <- p.adjust(MLE$p_value_purif_pos, method="fdr")
          
          # Permut Test 
          permut <- read.table(permut.file, h=T, row.names = 1)
          permut <- permut[match(rownames(MLE), rownames(permut)),]
          MLE$permut.pval <- permut$pval.high
          MLE$permut.FDR <- p.adjust(permut$pval.high, method="fdr")
          
          # Annotations
          MLE$species <- sp
          if (samp=="Schmidt10"){samp="Schmidt12"}
          MLE$TF <- TF
          MLE$sample <- samp
          MLE$sampleSize <- nrow(MLE)
          
          rownames(MLE) <- MLE$ID
          all_MLE_list[[paste0(sp, "_", TF)]] <- MLE
        }
      }
    }
    saveRDS(all_MLE_list, file=AllMerged)
  }else{all_MLE_list <- readRDS(AllMerged)
}
################################################################################
all_MLE <- do.call(rbind, all_MLE_list)

################################################################################
# Determine Conclusion according to FDR and alpha Tresholds
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_threshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_threshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
all_MLE$Conclusion <- factor(all_MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))

all_MLE$signif <- as.factor(ifelse(all_MLE$permut.pval <= pval_threshold, "Directional (+)", ifelse(all_MLE$permut.pval >= 1-pval_threshold, "Directional (-)", "Neutral")))
all_MLE$signif <- factor(all_MLE$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR <- "Neutral"
all_MLE$signif.FDR <- as.factor(ifelse(all_MLE$permut.FDR > fdr_threshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR <- factor(all_MLE$signif.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

################################################################################
all_MLE$length <- as.numeric(str_split_i(all_MLE$ID, ":", 3))-as.numeric(str_split_i(all_MLE$ID, ":", 2))
all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$length
all_MLE$species <- factor(all_MLE$species, levels=c("cat", "dog", "macaca", "human", "mouse",
                                                    "spretus", "rabbit", "caroli", "chicken", "rat", "drosophila"))


fwrite(all_MLE, paste0(path, "/../allMLE_vertebrates_", method, "_", minSub, "sub_all_exp.tsv"), sep="\t")
################################################################################