library(stringr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/"

species=c("caroli", "spretus", "dog", "cat", "human", "macaca", "rabbit", "rat", "mouse", "chicken")
sample=c("Wilson", "Myers", "Rensch", "Schmidt10", "Schmidt12", "Stefflova")
TFs=c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")
maxSub=150

all_MLE_list <- list()
for (sp in species){
  for (samp in sample){
    for (TF in TFs){
      MLE.file <- paste0(path, "MLE/", sp, "_", samp, "_", TF, "_MLE_summary_50bins.csv")
      delta.file <- paste0(path, "observed_deltas/", sp, "_", samp, "_", TF, "_deltas_ancestral_to_observed_deltaSVM.txt")
      permut.file <- paste0(path, "permutations_test/", sp, "_", samp, "_", TF, "_PosSelTest_deltaSVM_10000permutations.txt")
      
      if (file.exists(MLE.file)){
        print(paste(sp, samp, TF))
        
        # MLE Test
        MLE <- read.csv(MLE.file, h=T, row.names = 1)
        MLE$Conclusion <- as.factor(MLE$Conclusion)
        MLE$Conclusion <- factor(MLE$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
        MLE$ID <- rownames(MLE)
        MLE$species <- sp
        if (samp=="Schmidt10"){samp="Schmidt12"}
        MLE$TF <- TF
        MLE$sample <- samp
        
        # Observed deltas
        obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
        deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
        row.names(deltas) <- deltas$seq_name
        deltas <- deltas[row.names(MLE),]
        
        MLE$SVM <- deltas$SVM
        MLE$deltaSVM <- deltas$deltaSVM
        MLE$Diffdelta <- MLE$deltaSVM-MLE$SumObs
        
        # Permut Test 
        permut <- read.table(permut.file, h=T, row.names = 1)
        permut <- permut[rownames(MLE),]
        permut$signif <- as.factor(ifelse(permut$pval.high <= 0.01 | permut$pval.high >= 0.99, "Positive", "Neutral"))
        permut$signif <- factor(permut$signif, levels=c("Positive", "Neutral"))
        
        MLE <- cbind(MLE, permut[,c("pval.high", "signif")])
        all_MLE_list[[paste0(sp, "_", TF)]] <- MLE
      }
    }
  }
}

all_MLE <- do.call(rbind, all_MLE_list)
all_MLE$length <- as.numeric(str_split_i(all_MLE$ID, ":", 3))-as.numeric(str_split_i(all_MLE$ID, ":", 2))
all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$length

### Boxplots
all_MLE$species <- factor(all_MLE$species, levels=c("cat", "dog", "macaca", "human", "mouse",  "spretus", "rabbit", "caroli", "chicken", "rat"))

pdf(paste0(path, "figures/All_peaks_boxplots.pdf"), width = 10)
par(mfrow=c(1,2))
boxplot(all_MLE$Nmut~all_MLE$species, outline=F, xlab="", ylab="Substitutions", las=2)
boxplot(all_MLE$length~all_MLE$species, outline=F, xlab="", ylab="Length", las=2)

a <- boxplot(all_MLE$Prop.Sub~all_MLE$species, outline=F, xlab="", ylab="Substitution per base", las=2)
boxplot(all_MLE$Prop.Sub~all_MLE$TF, outline=F, xlab="TFs", ylab="Substitution per base", las=1)

boxplot(all_MLE$SumObs~all_MLE$species, outline=F, xlab="", ylab="Sum DeltaSVM Substitutions", las=2)
boxplot(all_MLE$SumObs~all_MLE$TF, outline=F, xlab="TFs", ylab="Sum DeltaSVM Substitutions", las=1)

boxplot(all_MLE$deltaSVM~all_MLE$species, outline=F, xlab="", ylab="DeltaSVM", las=2)
boxplot(all_MLE$deltaSVM~all_MLE$TF, outline=F, xlab="TFs", ylab="DeltaSVM", las=1)
dev.off()

#### Bar plots
col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Neutral model")

pdf(paste0(path, "figures/All_peaks_models_comparisons_barplots.pdf"), width = 11)
par(mfrow=c(1,2))

# Overlap Tests
conclu.test.Max <- prop.table(table(all_MLE$signif, all_MLE$Conclusion))
barplot(conclu.test.Max, col=c("forestgreen", "firebrick"), las=1, ylab="Proportion", main="maxLL vs Permutation Test")

conclu.test <- prop.table(table(all_MLE$Conclusion, all_MLE$signif))
barplot(conclu.test, col=col, las=1, ylab="Proportion", main="Permutation vs maxLL Test")

# Conclusion * Species
Conclu.sp <- prop.table(table(all_MLE$Conclusion, all_MLE$species), margin = 2)
barplot(Conclu.sp[1:2,], col=col, las=2, ylab="Proportion", main="MaxLL Test")

Conclu.sp.permut <- prop.table(table(all_MLE$signif, all_MLE$species), margin = 2)
barplot(Conclu.sp.permut[1,], col="forestgreen", las=2, ylab="Proportion", main="Permutation Test")
legend("topleft", legend=c("Positive", "Stabilising"), fill=col, bty="n")

# Conclusion * TFs
Conclu.TF <- prop.table(table(all_MLE$Conclusion, all_MLE$TF), margin = 2)
barplot(Conclu.TF[1:2,], col=col, las=1, ylab="Proportion", main="MaxLL Test")

Conclu.TF.permut <- prop.table(table(all_MLE$signif, all_MLE$TF), margin = 2)
barplot(Conclu.TF.permut[1,], col="forestgreen", las=1, ylab="Proportion", main="Permutation Test")
legend("topleft", legend=c("Positive", "Stabilising"), fill=col, bty="n")

# Conclusion * Sample
Conclu.samp <- prop.table(table(all_MLE$Conclusion, all_MLE$sample), margin = 2)
barplot(Conclu.samp[1:2,], col=col, las=1, ylab="Proportion", main="MaxLL Test")
legend("topright", legend=c("Positive", "Stabilising"), fill=col, bty="n")

Conclu.samp.permut <- prop.table(table(all_MLE$signif, all_MLE$sample), margin = 2)
barplot(Conclu.samp.permut[1,], col="forestgreen", las=1, ylab="Proportion", main="Permutation Test")

# Conclusion * Substitutions
all_MLE$Sub_class <- cut(all_MLE$Prop.Sub,  breaks=quantile(all_MLE$Prop.Sub, probs = seq(0, 1, 0.05)), include.lowest = T)
Conclu.Sub <- prop.table(table(all_MLE$Conclusion, all_MLE$Sub_class), margin = 2)
barplot(Conclu.Sub[1:2,], col=col, las=1, ylab="Proportion", xlab="Quantile Substitution per base",  main="MaxLL Test")

Conclu.Sub.permut <- prop.table(table(all_MLE$signif, all_MLE$Sub_class), margin = 2)
barplot(Conclu.Sub.permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab="Quantile Substitution per base",  main="Permutation Test")
legend("topleft", legend=c("Positive", "Stabilising"), fill=col, bty="n")

# Conclusion * SumObsSVM
all_MLE$SumObs_class <- cut(all_MLE$SumObs,  breaks=quantile(all_MLE$SumObs, probs = seq(0, 1, 0.05)), include.lowest = T)
Conclu.SumObs <- prop.table(table(all_MLE$Conclusion, all_MLE$SumObs_class), margin = 2)
barplot(Conclu.SumObs[1:2,], col=col, las=1, ylab="Proportion", xlab="Quantile Sum Observed DeltaSVM", main="MaxLL Test")

Conclu.SumObs.permut <- prop.table(table(all_MLE$signif, all_MLE$SumObs_class), margin = 2)
barplot(Conclu.SumObs.permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab="Quantile Sum Observed DeltaSVM",  main="Permutation Test")
legend("topleft", legend=c("Positive", "Stabilising"), fill=col, bty="n")

# Conclusion * deltaSVM
all_MLE$deltaSVM_class <- cut(all_MLE$deltaSVM,  breaks=quantile(all_MLE$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)
Conclu.deltaSVM <- prop.table(table(all_MLE$Conclusion, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:2,], col=col, las=1, ylab="Proportion", xlab="Quantile DeltaSVM", main="MaxLL Test")

Conclu.deltaSVM.permut <- prop.table(table(all_MLE$signif, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM.permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab="Quantile DeltaSVM", main="Permutation Test")


# Epistatic effect
boxplot(all_MLE$Diffdelta~all_MLE$Sub_class, outline=F, xlab="Quantile Substitution per base", ylab="Epistatic effect")

all_MLE$Diffdelta_class <- cut(all_MLE$Diffdelta,  breaks=quantile(all_MLE$Diffdelta, probs = seq(0, 1, 0.25)), include.lowest = T)
Conclu.epistatic <- prop.table(table(all_MLE$Conclusion, all_MLE$Diffdelta_class), margin = 2)
barplot(Conclu.epistatic[1:2,], col=col, las=2, ylab="Proportion", main="MaxLL Test")
legend("topright", legend=c("Positive", "Stabilising"), fill=col, bty="n")

dev.off()

### Heatmaps
# Substitution to deltaSVM
col <- colorRampPalette(brewer.pal(8, "Spectral"))(25)

all_MLE$Pos.bin <- ifelse(all_MLE$Conclusion=="Positive model", 1, 0)
all_MLE$Stab.bin <- ifelse(all_MLE$Conclusion=="Stabilizing model", 1, 0)
all_MLE$Rand.bin <- ifelse(all_MLE$Conclusion=="Neutral model", 1, 0)


# Heatmap 
# Calculate proportions
proportion <- as.data.frame(scale(t(table(all_MLE$SumObs_class, all_MLE$Sub_class))))
colnames(proportion) <- c("SumObs_class", "Sub_class", "Frequency")
proportion_stab <- all_MLE %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Stab.bin))
proportion_pos <- all_MLE %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Pos.bin))
proportion_rand <- all_MLE %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Rand.bin))

pdf(paste0(path, "figures/All_peaks_proportion_delta_sub_heatmap.pdf"), width = 10)
par(mfrow=c(1,1))
ggplot(proportion, aes(x = Sub_class, y = SumObs_class, fill = Frequency)) +
  geom_tile() +   scale_fill_gradientn(colors=rev(col)) + theme(legend.position = "none") + 
  labs(title = "Scaled frequency of peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_pos, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "forestgreen") + 
  labs(title = "Proportion of Positive peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_stab, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "deepskyblue3") + 
  labs(title = "Proportion of Stabilising peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_rand, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "firebrick") + 
  labs(title = "Proportion of Neutral peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")
dev.off()