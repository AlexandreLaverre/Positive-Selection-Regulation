# Analyse peaks droso
library(stringr)
library(qvalue)
library(progress)
library(data.table)
library(ROCR)


path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathPoly="/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/polymorphism_analyses/NarrowPeaks/drosophila/modERN/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
sp = "drosophila"
maxSub=150
minSub=2
fdr_threshold=0.05
pval_threshold=0.01

method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, "_", minSub, "sub_all_exp.Rds")
col=c("forestgreen", "orange", "deepskyblue4", "firebrick")
names(col) = c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")

################################################################################
# Get all data

if (file.exists(AllMerged)){
  # Load data
  all_MLE_list <- readRDS(AllMerged)
  
}else{
  # Get file paths
  exp.folders <- list.dirs(path, recursive = FALSE)
  exp_info <- read.csv("/Users/alaverre/Documents/Detecting_positive_selection/cluster/data/drosophila_modERN.csv", h=T)
  
  base_names <- paste(exp_info$Gene, exp_info$Strain, exp_info$Stage, sep = "_")
  name_counts <- ave(base_names, base_names, FUN = seq_along)
  unique_names <- paste(base_names, name_counts, sep = "_")
  exp_info$exp <- unique_names
  
  pb <- progress_bar$new(format = "  Processing [:bar] :percent (:elapsed s) | ETA :eta", total = length(exp.folders), clear = FALSE, width = 60)
  
  all_MLE_list <- list()
  
  for (i in 1:length(exp.folders)){
    exp.path = exp.folders[i]
    exp <- basename(exp.path)
    pb$tick()
    
    message(exp)
    if (exp %in% names(all_MLE_list)){next}
    if (exp=="log"){next}
    
    # Data
    delta.file <- paste0(exp.path, "/deltas/ancestral_to_observed_deltaSVM.txt")
    MLE.file <- paste0(exp.path, "/Tests/MLE_summary_", method, ".csv")
    permut.file <- paste0(exp.path, "/Tests/PosSelTest_deltaSVM_10000permutations_two_tailed_ancestral.txt")
    stats.file <- paste0(exp.path, "/sequences/focal_ancestral_substitutions_stats.txt")
    polymorphism.file <- paste0(pathPoly, exp, "/overlap_peaks.txt")
    
    # MaxLL Test and FDR correction
    MLE <- read.csv(MLE.file, h=T, row.names = 1)
    MLE$exp <- exp
    MLE$DevStage <-  sub(".*_(.*_.*)$", "\\1", exp)
    MLE$TF <- sub("_.*", "", exp)
    
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
    MLE$Diffdelta <- MLE$deltaSVM-MLE$SumObs
    
    # Filter out peaks with >20% indel
    MLE$ID <- gsub("[_\\:]([[:alpha:]]).*_.*", "", rownames(MLE))
    MLE$original_length <- as.numeric(str_split_i(MLE$ID, ":", 3))-as.numeric(str_split_i(MLE$ID, ":", 2))
    MLE$indel <- MLE$original_length - seq_stats$Length
    MLE$Prop.indel <- (MLE$indel*100)/MLE$original_length
    MLE$original_sample_size <- nrow(MLE)
    MLE <- MLE[which(MLE$Prop.indel<20),]
    
    # Substitution filter
    MLE <- MLE[which(MLE$Nmut >= minSub),]
    
    if (nrow(MLE)<1){next}
    MLE$filtered_sample_size <- nrow(MLE)
    seq_stats <- seq_stats[row.names(MLE),]
    
    # FDR correction
    MLE$FDR_null_pos <- p.adjust(MLE$p_value_null_pos, method="fdr")
    MLE$FDR_null_purif <- p.adjust(MLE$p_value_null_purif, method="fdr")
    MLE$FDR_purif_pos <- p.adjust(MLE$p_value_purif_pos, method="fdr")
    
    # Polymorphism data
    polymorphism <- read.table(polymorphism.file)
    polymorphism <- unique(polymorphism[,1:2])
    NbSNP <- as.data.frame(table(polymorphism$V1))
    rownames(NbSNP) <- NbSNP$Var1
    if (grepl("^chr", rownames(NbSNP)[1])) {rownames(NbSNP) <- sub("^chr[^_]*_", "", rownames(NbSNP))}
    
    MLE$NbSNP <- NbSNP[rownames(MLE),]$Freq
    MLE$NbSNP[is.na(MLE$NbSNP)] <- 0 # When SNP not found in polymorphism file
    if (sum(MLE$NbSNP)==0){message("No SNP found in polymorphism file for ", exp)}
    
    # Permutations Test
    permut <- read.table(permut.file, h=T, row.names = 1)
    permut <- permut[rownames(MLE),]
    permut$FDR <- p.adjust(permut$pval.high, method="fdr")
    permut$FDR.two <- p.adjust(permut$pval.two.tailed, method="fdr")
    
    MLE <- cbind(MLE, permut[,c("pval.high", "FDR", "FDR.two")], seq_stats[,c("Length", "GC_Content", "Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak")])
    all_MLE_list[[exp]] <- MLE
    
  }

  # Save output 
  saveRDS(c(all_MLE_list), file=AllMerged)
}

################################################################################
all_MLE <- do.call(rbind, all_MLE_list)
all_MLE$peakID <- unlist(lapply(all_MLE_list, rownames))
rownames(all_MLE) <- all_MLE$peakID

# Filter Low Sample Size
all_MLE2 <- all_MLE[which(all_MLE$original_sample_size > 1500 & all_MLE$AUC >0.8),]
all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$Length

################################################################################
# Determine Conclusion according to FDR and alpha Tresholds
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_threshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_threshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
all_MLE$Conclusion <- factor(all_MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))

all_MLE$signif <- as.factor(ifelse(all_MLE$pval.high <= pval_threshold, "Directional (+)", ifelse(all_MLE$pval.high >= 1-pval_threshold, "Directional (-)", "Neutral")))
all_MLE$signif <- factor(all_MLE$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR <- "Neutral"
all_MLE$signif.FDR <- as.factor(ifelse(all_MLE$FDR > fdr_threshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR <- factor(all_MLE$signif.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR.two <- "Neutral"
all_MLE$signif.FDR.two <- as.factor(ifelse(all_MLE$FDR.two > fdr_threshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR.two <- factor(all_MLE$signif.FDR.two, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

fwrite(all_MLE, paste0(path, "/../allMLE_drosophila_", method, "_", minSub, "sub_all_exp.tsv"), sep="\t")

################################################################################