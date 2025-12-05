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
  exp_info <- read.csv("/Users/alaverre/Documents/Detecting_positive_selection/cluster/data/droso_exp_info.csv", h=T)
  
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
    
    # Quality filtered 
    #MLE <- MLE[which(MLE$SVM > -350),]
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
all_MLE <- all_MLE[which(all_MLE$original_sample_size > 1500 & all_MLE$AUC >0.8),]
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
all_MLE$Sub_class <- cut(all_MLE$Nmut,  breaks=c(0,2,5,10,15,20,25,50,100,150), include.lowest = T,
                         labels = c("2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-50", "51-100", ">100"))

all_MLE$Prop_Sub_class <- cut(all_MLE$Prop.Sub,  breaks=quantile(all_MLE$Prop.Sub, probs = seq(0, 1, 0.1)), include.lowest = T, dig.lab=1)
all_MLE$Prop_indel_class <- cut(all_MLE$Prop.indel,  breaks=seq(0, 100, 10), include.lowest = T,
                                labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"))

################################################################################
par(mfrow=c(2,2), mar=c(4,4,3,1))
hist(all_MLE$Nmut, breaks=1000, xlab="Nb Substitution", main="")
hist(all_MLE$Length, breaks=100, xlab="Peak Length (no indel)", main="")
boxplot(all_MLE$Length~all_MLE$Sub_class, outline=F, xlab="Nb Substitution", ylab="Peaks length (no indel)", notch=T, main="")
hist(all_MLE$Prop.Sub, breaks=200, xlab="Substitution per base", main="")
boxplot(all_MLE$Length~all_MLE$Prop_Sub_class, outline=F, xlab="Substitution per base", ylab="Peaks length (no indel)", notch=T, main="")

# Indel
par(mfrow=c(2,2), mar=c(4,4,3,1))
hist(all_MLE$Prop.indel, breaks=100, xlab="% indel", main="")
boxplot(all_MLE$Prop.Sub~all_MLE$Prop_indel_class, outline=F, ylab="Substitution per base", xlab="% indel", notch=T, main="")
boxplot(all_MLE$original_length~all_MLE$Prop_indel_class, outline=F, xlab="% indel", ylab="Original Length", notch=T, main="")
boxplot(all_MLE$Length~all_MLE$Prop_indel_class, outline=F, xlab="% indel", ylab="Length no indel", notch=T, main="")

all_MLE$DevStage <- gsub("_.*", "", all_MLE$DevStage)
all_MLE[grepl("instarlarva", all_MLE$DevStage),]$DevStage <- "larva"
all_MLE$DevStage <- factor(all_MLE$DevStage, levels=c("embryonic", "larva", "prepupa", "pupa", "adult", "unknown"))

################################################################################
# Correlation with features
#pdf(paste0(pathFigures, "drosophila_peaks_tmp.pdf"), width=10, height=5)
all_MLE$filtered_sample_size <- as.factor(all_MLE$filtered_sample_size)
par(mfrow=c(1,1))
# Sample Size effect on SVM
meanSVM <- tapply(all_MLE$SVM, all_MLE$filtered_sample_size, median)
plot(meanSVM, xlab="SampleSize", ylab="median SVM", type="b", main="Drosophila peaks SVM scores", pch=19, cex=0.5, xaxt = "n")
index = round(seq(1, length(meanSVM), length.out = 10))
axis(1, at = index, labels = names(meanSVM)[index], las = 1) 
min_samp=length(which(as.numeric(levels(all_MLE$original_sample_size))<2500))
abline(v=min_samp, col="red", lwd=1.5)
abline(h=-350, col="red", lwd=1.5)
mtext(paste("Add threshold \n N>2500 peaks?"), col="red", side=1, line=-2, cex=0.9)

################################################################################
# Sample Size on proportion
par(mfrow=c(2,2), mar=c(4,4,3,1))
# maxLL - Sample Size

Conclu.sample <- prop.table(table(all_MLE$Conclusion, all_MLE$filtered_sample_size) , margin = 2)
bar_pos <- barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="maxLL Test (pval< 1%)", xaxt = "n")
min_samp=length(which(as.numeric(levels(all_MLE$filtered_sample_size))<2500))
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

Conclu.sample.FDR <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$filtered_sample_size) , margin = 2)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="maxLL Test (FDR<5%)", xaxt = "n")
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

# Permut - Sample Size
Conclu.sample <- prop.table(table(all_MLE$signif, all_MLE$filtered_sample_size) , margin = 2)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="Permut (pval< 1%)", xaxt = "n")
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

Conclu.sample.FDR <- prop.table(table(all_MLE$signif.FDR, all_MLE$filtered_sample_size) , margin = 2)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="Permut (FDR<10%)", xaxt = "n")
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

################################################################################
# Sample Size on total number
par(mfrow=c(2,2), mar=c(4,4,3,1))
# maxLL - Sample Size
all_MLE$filtered_sample_size <- as.factor(all_MLE$filtered_sample_size)
Conclu.sample <- table(all_MLE$Conclusion, all_MLE$filtered_sample_size)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="maxLL Test (pval< 1%)")
min_samp=length(which(as.numeric(levels(all_MLE$filtered_sample_size))<2500))
abline(v=min_samp, col="red")

Conclu.sample.FDR <-table(all_MLE$FDR_Conclusion, all_MLE$filtered_sample_size)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="maxLL Test (FDR<10%)")
abline(v=min_samp, col="red")

# Permut - Sample Size
Conclu.sample <- table(all_MLE$signif, all_MLE$filtered_sample_size)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="Permut (pval< 1%)")
min_samp=length(which(as.numeric(levels(all_MLE$filtered_sample_size))<2500))
abline(v=min_samp, col="red")

Conclu.sample.FDR <- table(all_MLE$signif.FDR, all_MLE$filtered_sample_size)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="Permut (FDR<10%)")
abline(v=min_samp, col="red")

################################################################################
# Delta SVM - maxLL
all_MLE$deltaSVM_class <- cut(all_MLE$deltaSVM,  breaks=quantile(all_MLE$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)

Conclu.deltaSVM <- prop.table(table(all_MLE$Conclusion, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:3,], col=col, las=1, ylab="Proportion Directional", xlab=expression(Delta*"SVM (quantile)"), main="maxLL (pval<1%)")
legend("topleft", legend=c("Directional (+)", "Directional (-)", "Stabilising"), fill=col, bty="n")

Conclu.deltaSVM <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:3,], col=col, las=1, ylab="Proportion Directional", xlab="Phenotypic Change (quantile)", main="maxLL (FDR < 0.1)")

# Delta SVM - Permut
Conclu.deltaSVM <- prop.table(table(all_MLE$signif, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:3,], col=col, las=1, ylab="Proportion Directional", xlab=expression(Delta*"SVM (quantile)"), main="Permut (pval<1%)")

Conclu.deltaSVM <- prop.table(table(all_MLE$signif.FDR, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:3,], col=col, las=1, ylab="Proportion Directional", xlab=expression(Delta*"SVM (quantile)"), main="Permut (FDR<10%)")

################################################################################
# Developmental Stage
Conclu.dev.FDR <- prop.table(table(all_MLE$FDR_Conclusion, as.factor(all_MLE$DevStage)), margin = 2)
barplot(Conclu.dev.FDR[1:2,1:5], col=col, las=2, ylab="Proportion Directional", main="maxLL (FDR < 0.1)")

Conclu.dev.FDR <- prop.table(table(all_MLE$signif.FDR, as.factor(all_MLE$DevStage)), margin = 2)
barplot(Conclu.dev.FDR[1:2,1:5], col=col, las=2, ylab="Proportion Directional", main="Permut (FDR < 0.1)")

# Substitutions
all_MLE$Sub_class <- cut(all_MLE$Nmut,  breaks=c(0,2,5,10,15,20,25,50,100,150), include.lowest = T,
                         labels = c("2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-50", "51-100", ">100"))

Conclu.Sub.FDR <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$Sub_class), margin = 2)
barplot(Conclu.Sub.FDR[1:2,], col=col, las=2,  ylab="Proportion Directional", xlab="Nb Sub", main="maxLL (FDR < 0.1)")

Conclu.Sub.FDR <- prop.table(table(all_MLE$signif.FDR, all_MLE$Sub_class), margin = 2)
barplot(Conclu.Sub.FDR[1:2,], col=col, las=2,  ylab="Proportion Directional", xlab="Nb Sub", main="Permut (FDR < 0.1)")


# GC content
quantile_5 = paste0(seq(5, 100, 5), "%")
all_MLE$GC_class <- cut(all_MLE$GC_Content,  breaks=quantile(all_MLE$GC_Content, probs = seq(0, 1, 0.05)), include.lowest = T)

Conclu.GC_class <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$GC_class), margin = 2)
barplot(Conclu.GC_class[1:2,], col=col, las=1, ylab="Proportion Directional", xlab="GC Content (quantile)", main="MaxLL (FDR)", names=quantile_5)

Conclu.GC_class_permut <- prop.table(table(all_MLE$signif.FDR, all_MLE$GC_class), margin = 2)
barplot(Conclu.GC_class_permut[1:2,], col=col, las=1, ylab="Proportion Directional", xlab="GC Content (quantile)", main="Permut (FDR)", names=quantile_5)
legend("topright", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")


par(mfrow=c(1,2))
for (sub_type in c("Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak")){
  print(sub_type)
  all_MLE[[paste0(sub_type, "_prop")]] <- all_MLE[[sub_type]]/all_MLE$Nmut
  sub_type = paste0(sub_type, "_prop")
  all_MLE$sub_class <- cut(all_MLE[[sub_type]],  breaks= seq(0,1,0.1), include.lowest = T)
  Conclu.sub_class <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$sub_class), margin = 2)
  barplot(Conclu.sub_class[1:2,], col=col, las=1, ylab="Proportion", xlab=paste("Quantile", sub_type), main="MaxLL")
  
  Conclu.sub_class_permut <- prop.table(table(all_MLE$signif.FDR, all_MLE$sub_class), margin = 2)
  barplot(Conclu.sub_class_permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab=paste("Quantile", sub_type), main="Permutation")
}

#dev.off()

################################################################################