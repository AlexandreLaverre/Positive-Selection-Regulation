# Analyse peaks droso
library(stringr)
library(qvalue)
library(progress)


path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathPoly="/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/polymorphism_analyses/NarrowPeaks/drosophila/modERN/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
sp = "drosophila"
maxSub=150
alpha_treshold=0.01

method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, ".Rds")
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
  
  base_names <- paste(exp_info$Gene, exp_info$Stage, sep = "_")
  name_counts <- ave(base_names, base_names, FUN = seq_along)
  unique_names <- paste(base_names, name_counts, sep = "_")
  exp_info$ID <- unique_names
  
  pb <- progress_bar$new(format = "  Processing [:bar] :percent (:elapsed s) | ETA :eta", total = length(exp.folders), clear = FALSE, width = 60)
  
  all_MLE_list <- list()
  for (i in 1:length(exp.folders)){
    exp.path = exp.folders[i]
    exp <- basename(exp.path)
    pb$tick()
    
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
    
    # Filter out peaks with >20% indel
    MLE$ID <- gsub("[_\\:]([[:alpha:]]).*_.*", "", rownames(MLE))
    MLE$original_length <- as.numeric(str_split_i(MLE$ID, ":", 3))-as.numeric(str_split_i(MLE$ID, ":", 2))
    MLE$indel <- MLE$original_length - seq_stats$Length
    MLE$Prop.indel <- (MLE$indel*100)/MLE$original_length
    MLE <- MLE[which(MLE$Prop.indel<20),]
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
    MLE$NbSNP <- NbSNP[rownames(MLE),]$Freq
    MLE$NbSNP[is.na(MLE$NbSNP)] <- 0 # When SNP not found in polymorphism file
    
    # Permutations Test
    permut <- read.table(permut.file, h=T, row.names = 1)
    permut <- permut[rownames(MLE),]
    permut$FDR <- p.adjust(permut$pval.high, method="fdr")
    permut$FDR.two <- p.adjust(permut$pval.two.tailed, method="fdr")
    
    # Observed deltas
    obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
    deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
    row.names(deltas) <- deltas$seq_name
    deltas <- deltas[row.names(MLE),]
    
    MLE$SVM <- deltas$SVM
    MLE$deltaSVM <- deltas$deltaSVM
    MLE$Diffdelta <- MLE$deltaSVM-MLE$SumObs
    
    MLE <- cbind(MLE, permut[,c("pval.high", "FDR", "FDR.two")], seq_stats[,c("Length", "GC_Content", "Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak")])
    all_MLE_list[[exp]] <- MLE
    
    # Model performance and Sample Size
    cv<-fread(paste0(exp.path, "/Model/", exp, ".cvpred.txt"))
    colnames(cv)<-c("position","prediction","real_state","cv_number")
    pred <- prediction(cv$prediction, cv$real_state) 
    auc_result <- performance(pred, measure = "auc")
    AUC = signif(unlist(slot(auc_result, "y.values")), digits=3)
    
    all_AUC$exp[i] <- exp
    all_AUC$AUC[i] <- AUC
    all_AUC$SampleSize[i] <- nrow(MLE)
    
  }
      
  # Save output 
  saveRDS(c(all_MLE_list, all_AUC), file=AllMerged)
}

################################################################################

all_MLE <- do.call(rbind, all_MLE_list)
all_MLE$peakID <- unlist(lapply(all_MLE_list, rownames))
rownames(all_MLE) <- all_MLE$peakID

# Filter Low Sample Size
all_MLE <- all_MLE[which(all_MLE$SampleSize > 2500),]

all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$Length
all_MLE$Sub_class <- cut(all_MLE$Nmut,  breaks=c(0,2,5,10,15,20,25,50,100,150), include.lowest = T,
                         labels = c("2", "3-5", "6-10", "11-15", "16-20", "21-25", "26-50", "51-100", ">100"))

all_MLE$Prop_Sub_class <- cut(all_MLE$Prop.Sub,  breaks=quantile(all_MLE$Prop.Sub, probs = seq(0, 1, 0.1)), include.lowest = T, dig.lab=1)
all_MLE$Prop_indel_class <- cut(all_MLE$Prop.indel,  breaks=seq(0, 100, 10), include.lowest = T,
                                labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"))


par(mfrow=c(2,2), mar=c(4,4,3,1))
#hist(all_MLE$Nmut, breaks=1000, xlab="Nb Substitution", main="")
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
# Determine Conclusion according to FDR and alpha Tresholds
fdr_treshold=0.1
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_treshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_treshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
all_MLE$Conclusion <- factor(all_MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))

all_MLE$signif <- as.factor(ifelse(all_MLE$pval.high <= alpha_treshold, "Directional (+)", ifelse(all_MLE$pval.high >= 0.99, "Directional (-)", "Neutral")))
all_MLE$signif <- factor(all_MLE$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR <- "Neutral"
all_MLE$signif.FDR <- as.factor(ifelse(all_MLE$FDR > fdr_treshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR <- factor(all_MLE$signif.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR.two <- "Neutral"
all_MLE$signif.FDR.two <- as.factor(ifelse(all_MLE$FDR.two > fdr_treshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR.two <- factor(all_MLE$signif.FDR.two, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

# Write table
fwrite(all_MLE, paste0(path, "/../allMLE_drosophila_", method, "_filtered.tsv"), sep="\t")

################################################################################
# Correlation with features
#pdf(paste0(pathFigures, "drosophila_peaks_tmp.pdf"), width=10, height=5)
all_MLE$SampleSize <- as.factor(all_MLE$SampleSize)
par(mfrow=c(1,1))
# Sample Size effect on SVM
meanSVM <- tapply(all_MLE$SVM, all_MLE$SampleSize, median)
plot(meanSVM, xlab="SampleSize", ylab="median SVM", type="b", main="Drosophila peaks SVM scores", pch=19, cex=0.5, xaxt = "n")
index = round(seq(1, length(meanSVM), length.out = 10))
axis(1, at = index, labels = names(meanSVM)[index], las = 1) 
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red", lwd=1.5)
mtext(paste("Add threshold \n N>2500 peaks?"), col="red", side=1, line=-2, cex=0.9)

################################################################################
# Sample Size on proportion
par(mfrow=c(2,2), mar=c(4,4,3,1))
# maxLL - Sample Size

Conclu.sample <- table(all_MLE$Conclusion, all_MLE$SampleSize)
bar_pos <- barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="maxLL Test (pval< 1%)", xaxt = "n")
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

Conclu.sample.FDR <- table(all_MLE$FDR_Conclusion, all_MLE$SampleSize)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="maxLL Test (FDR<10%)", xaxt = "n")
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

# Permut - Sample Size
Conclu.sample <- prop.table(table(all_MLE$signif, all_MLE$SampleSize), margin = 2)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="Permut (pval< 1%)", xaxt = "n")
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

Conclu.sample.FDR <- prop.table(table(all_MLE$signif.FDR, all_MLE$SampleSize), margin = 2)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Proportion Directional", xlab="Sample Size", main="Permut (FDR<10%)", xaxt = "n")
abline(v=min_samp, col="red")
axis(1, at = bar_pos[index], labels = names(meanSVM)[index], las = 2) 

################################################################################
# Sample Size on total number
par(mfrow=c(2,2), mar=c(4,4,3,1))
# maxLL - Sample Size
all_MLE$SampleSize <- as.factor(all_MLE$SampleSize)
Conclu.sample <- table(all_MLE$Conclusion, all_MLE$SampleSize)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="maxLL Test (pval< 1%)")
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red")

Conclu.sample.FDR <-table(all_MLE$FDR_Conclusion, all_MLE$SampleSize)
barplot(Conclu.sample.FDR[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="maxLL Test (FDR<10%)")
abline(v=min_samp, col="red")

# Permut - Sample Size
Conclu.sample <- table(all_MLE$signif, all_MLE$SampleSize)
barplot(Conclu.sample[1:2,], col=col, las=2, ylab="Nb peaks Directional", xlab="Sample Size", main="Permut (pval< 1%)")
min_samp=length(which(as.numeric(levels(all_MLE$SampleSize))<2500))
abline(v=min_samp, col="red")

Conclu.sample.FDR <- table(all_MLE$signif.FDR, all_MLE$SampleSize)
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
## Polymorphism
TFS = levels(as.factor(all_MLE$exp))
fisher_df <- data.frame(TF = TFS, odds_ratio = NA, p_value = NA)

for (i in seq_along(TFS)){
  TF <- TFS[i]
  print(TF)
  sub_MLE <- all_MLE[which(all_MLE$exp==TF),]
  # Directional vs Other
  pos <- sub_MLE[grepl("Directional", sub_MLE$FDR_Conclusion),]
  nonPos <- sub_MLE[!grepl("Directional", sub_MLE$FDR_Conclusion),]
  subNumb <-c(sum(pos$Nmut),sum(nonPos$Nmut))
  polyNumb <-c(sum(pos$NbSNP),sum(nonPos$NbSNP))
  
  test <- fisher.test(matrix(c(subNumb,polyNumb),nrow = 2,ncol = 2))
  
  fisher_df$odds_ratio[i] <- test$estimate
  fisher_df$p_value[i] <- test$p.value
}

# Volcano Plot
par(mfrow=c(1,1))
fisher_df <- fisher_df[which(fisher_df$odds_ratio>0),]
rownames(fisher_df) <- fisher_df$TF
fisher_df$SampleSize <- all_MLE$SampleSize[match(fisher_df$TF, all_MLE$exp)]
plot(fisher_df$odds_ratio, -log10(fisher_df$p_value),
     pch=20, cex = 0.8, col=ifelse(fisher_df$p_value < 0.05, ifelse(fisher_df$odds_ratio>1, "red", "deepskyblue3"), "grey"),
     xlab = "Odds Ratio", ylab = "-log10(p-value)", cex.lab=1.4, cex.lab=1.2,
     main = "maxLL (FDR)")
legend("topleft", legend=c("Directional > Other", "Other > Directional", "NS"), col=c("red", "deepskyblue3", "grey"), bty="n",  pch=20)


# Map values to colors
cols <- c("#440154", "#3b528b", "#21908d", "#5dc963", "#fde725")
val_scaled <- as.numeric(cut(fisher_df$SampleSize, breaks=c(0,3000,5000,10000,15000,max(fisher_df$SampleSize))))  # scale values into 100 bins
fisher_df$point_colors <- cols[val_scaled]

plot(fisher_df$odds_ratio, -log10(fisher_df$p_value),
     pch=20, cex = 0.8, col=fisher_df$point_colors,
     xlab = "Odds Ratio", ylab = "-log10(p-value)", cex.lab=1.4, cex.lab=1.2,
     main = "maxLL (FDR)")

legend_labels <- c("2.5–3k", "3k–5k", "5k–10k", "10k–15k", "15k+")

# Add the legend
legend("topright", legend = legend_labels, col = cols, title = "Sample Size", border = NA, cex = 0.8,pt.cex = 1.5, bty="n", pch = 20)


# Simple Barplot
par(mar=c(4,5,3,1))
ymax=1.6
bp<-barplot(subNumb/polyNumb, ylim=c(0,ymax), main="D1_pupa", cex.lab=1.4, cex.main=1.5, cex.axis =1.2, cex.names=1.5,
            names=c("Directional","Other"), ylab="# Substitutions / # polymorphisms", col=c("red", "deepskyblue3"))

text(x=bp, y=0+1.5/15, labels=paste0("n=", c(nrow(pos),nrow(nonPos))),cex=1.5)
fTest<-fisher.test(matrix(c(subNumb,polyNumb),nrow = 2,ncol = 2))
legend("topright",legend=paste0("Ratio=", signif(fTest$estimate,2), "; p=",signif(fTest$p.value,2)),bty = 'n',cex=1.5)

################################################################################
##### AUC 
all_AUC <- all_AUC[!is.na(all_AUC$AUC),]
rownames(all_AUC) <- all_AUC$exp

# Plots
hist(all_AUC$AUC, breaks=100, xlab="AUC", main="Drosophila - all datasets")
plot(all_AUC$AUC~all_AUC$SampleSize, cex=0.8, xlab="Sample Size", ylab="AUC", pch=19)
abline(v=1500, col="red")

################################################################################
