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
alpha_treshold=0.01
qval_treshold=0.2

col=c("forestgreen", "orange", "deepskyblue4", "firebrick")
names(col) = c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")

col_3=c("forestgreen", "deepskyblue4", "firebrick")
names(col_3) = c("Directional", "Stabilizing", "Neutral model")

method="exact" # quantile_50bins_threshold_0.01
################################################################################
# Get an object of all merged results
AllMerged = paste0(path, "allMLE_list_", method, ".Rds")
if (!file.exists(AllMerged)){
  all_MLE_list <- list()
  for (sp in species){
    for (samp in sample){
      for (TF in TFs){
        MLE.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/MLE_summary_", method, ".csv")
        delta.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/deltas/ancestral_to_observed_deltaSVM.txt")
        permut.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
        stats.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/sequences/focal_substitutions_stats.txt")
        
        if (file.exists(MLE.file)){
          print(paste(sp, samp, TF))
          
          # MLE Test
          MLE <- read.csv(MLE.file, h=T, row.names = 1)
          #MLE <- MLE[which(MLE$Nmut>4),]
          if (nrow(MLE)<500){
            print(paste(sp, samp, TF, "not enough sequences."))
            next}
          MLE$Conclusion <- ifelse(MLE$p_value_null_pos<alpha_treshold & MLE$p_value_null_purif>alpha_treshold & MLE$p_value_purif_pos>alpha_treshold, 
                                   ifelse(MLE$AlphaPos>MLE$BetaPos,"Directional (+)", "Directional (-)"), MLE$Conclusion)
          MLE$Conclusion <- ifelse(MLE$Conclusion == "Neutral model", "Neutral", MLE$Conclusion)
          MLE$Conclusion <- as.factor(MLE$Conclusion)
          MLE$Conclusion <- factor(MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          
          MLE$FDR_null_pos <- p.adjust(MLE$p_value_null_pos, method="fdr")
          MLE$FDR_null_purif <- p.adjust(MLE$p_value_null_purif, method="fdr")
          MLE$FDR_purif_pos <- p.adjust(MLE$p_value_purif_pos, method="fdr")
          MLE$FDR_Conclusion <- "Neutral"
          MLE$FDR_Conclusion <- ifelse(MLE$FDR_null_purif <= qval_treshold, "Stabilizing", MLE$FDR_Conclusion)
          MLE$FDR_Conclusion <- ifelse(MLE$FDR_purif_pos <= qval_treshold, ifelse(MLE$AlphaPos>MLE$BetaPos,"Directional (+)", "Directional (-)"), MLE$FDR_Conclusion)
          MLE$FDR_Conclusion <- factor(MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          
          MLE$qval_null_pos <- qvalue(MLE$p_value_null_pos)$qvalues
          MLE$qval_null_purif <- qvalue(MLE$p_value_null_purif)$qvalues
          MLE$qval_purif_pos <- qvalue(MLE$p_value_purif_pos)$qvalues
          MLE$qval_Conclusion <- "Neutral"
          MLE$qval_Conclusion <- ifelse(MLE$qval_null_purif <= qval_treshold, "Stabilizing", MLE$qval_Conclusion)
          MLE$qval_Conclusion <- ifelse(MLE$qval_purif_pos <= qval_treshold, ifelse(MLE$AlphaPos>MLE$BetaPos,"Directional (+)", "Directional (-)"), MLE$qval_Conclusion)
          MLE$qval_Conclusion <- factor(MLE$qval_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          
          MLE$species <- sp
          if (samp=="Schmidt10"){samp="Schmidt12"}
          MLE$TF <- TF
          MLE$sample <- samp
          MLE$sampleSize <- nrow(MLE)
          
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
          
          # Permut Test 
          permut <- read.table(permut.file, h=T, row.names = 1)
          permut <- permut[rownames(MLE),]
          permut$FDR <- p.adjust(permut$pval.high, method="fdr")
          permut$signif.FDR <- "Neutral"
          permut$signif.FDR <- as.factor(ifelse(permut$FDR > qval_treshold, "Neutral", ifelse(permut$deltaSVM > 0, "Directional (+)", "Directional (-)")))
          permut$signif.FDR <- factor(permut$signif.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          permut$signif <- as.factor(ifelse(permut$pval.high <= alpha_treshold, "Directional (+)", ifelse(permut$pval.high >= 0.99, "Directional (-)", "Neutral")))
          permut$signif <- factor(permut$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

          # Permut two-tailed
          permut$pval.two.tailed <- apply(permut, 1, function(x) 2*min(as.numeric(x["pval.high"]), 1-as.numeric(x["pval.high"])))
          permut$FDR.two <- p.adjust(permut$pval.two.tailed, method="fdr")
          #permut$qval.two <- qvalue(permut$pval.two.tailed)$qvalues
          
          permut$signif.two <- as.factor(ifelse(permut$pval.two.tailed > alpha_treshold, "Neutral", ifelse(permut$deltaSVM > 0, "Directional (+)", "Directional (-)")))
          permut$signif.two <- factor(permut$signif.two, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          permut$signif.two.FDR <- "Neutral"
          permut$signif.two.FDR <- as.factor(ifelse(permut$FDR.two > qval_treshold, "Neutral", ifelse(permut$deltaSVM > 0, "Directional (+)", "Directional (-)")))
          permut$signif.two.FDR <- factor(permut$signif.two.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          #permut$signif.two.qval <- "Neutral"
          #permut$signif.two.qval <- as.factor(ifelse(permut$qval.two > qval_treshold, "Neutral", ifelse(permut$deltaSVM > 0, "Directional (+)", "Directional (-)")))
          #permut$signif.two.qval <- factor(permut$signif.two.qval, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
          
          #correct error in IDs
          rownames(MLE) <- gsub("_.*:\\d+:\\d+:", "_", rownames(MLE))
          MLE$peaks_ID <- gsub(".*:\\d+:\\d+[_\\:]", "", rownames(MLE))
          MLE$ID <- gsub("[_\\:]([[:alpha:]]).*_.*", "", rownames(MLE))
          
          MLE <- cbind(MLE, permut[,c("pval.high", "FDR", "pval.two.tailed", "signif", "signif.FDR", "signif.two.FDR")], seq_stats[,c("GC_Content", "Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak")])
          all_MLE_list[[paste0(sp, "_", TF)]] <- MLE
        }
      }
    }
  }
  saveRDS(all_MLE_list, file=AllMerged)
  
}else{
  all_MLE_list <- readRDS(AllMerged)
}
################################################################################

all_MLE <- do.call(rbind, all_MLE_list)

all_MLE$length <- as.numeric(str_split_i(all_MLE$ID, ":", 3))-as.numeric(str_split_i(all_MLE$ID, ":", 2))
all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$length
all_MLE$species <- factor(all_MLE$species, levels=c("cat", "dog", "macaca", "human", "mouse",
                                                    "spretus", "rabbit", "caroli", "chicken", "rat", "drosophila"))

species=species[species != "drosophila"]
all_MLE <- all_MLE[which(all_MLE$species != "drosophila"),]

# Remove unusual chromosomes 
all_MLE <- all_MLE[which(!is.na(all_MLE$pval.high)),]

# Correlation matrix
selected_vars <- c("SVM", "deltaSVM", "Nmut", "pval.high", "pval.two.tailed", "FDR",
                   "LL_pos", "LL_purif", "LL_neutral", "p_value_null_pos", "p_value_null_purif", "p_value_purif_pos")


cor_matrix <- cor(all_MLE[,selected_vars], use = "pairwise.complete.obs", method = "pearson")
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black")

num_vars <- all_MLE[sapply(all_MLE, is.numeric)]
cor_matrix <- cor(num_vars, use = "pairwise.complete.obs")
cor_matrix["Nmut", ]
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black")

library(effectsize)

eta2_results <- sapply(all_MLE[,selected_vars], function(x) {
  model <- aov(x ~ all_MLE$sampleSize)
  eta_squared(model)$Eta2[1]
})

################################################################################

all_MLE$diffTests <- ifelse(all_MLE$Conclusion=="Directional (+)" & all_MLE$signif.FDR=="Neutral", "Dir/Null", 
                            ifelse(all_MLE$Conclusion=="Neutral" & all_MLE$signif.FDR=="Directional (+)", "Null/Dir", 
                                   ifelse(all_MLE$Conclusion=="Neutral" & all_MLE$signif.FDR=="Neutral", "Null/Null", 
                                          ifelse(all_MLE$Conclusion=="Directional (+)" & all_MLE$signif.FDR=="Directional (+)", "Dir/Dir", "Other"))))

all_MLE$diffTests <- factor(all_MLE$diffTests, levels=c("Null/Null", "Dir/Dir", "Dir/Null", "Null/Dir", "Other"))
par(mfrow=c(2,2))
boxplot(all_MLE$Nmut~all_MLE$qval_Conclusion, outline=F, notch=T, xlab="MLL Conclusion", ylab="Substitutions", ylim=c(0,35))
boxplot(all_MLE$Nmut~all_MLE$signif.FDR, outline=F, notch=T, xlab="Permut Conclusion", ylab="Substitutions")

boxplot(all_MLE$Nmut~all_MLE$diffTests, outline=F, notch=T, xlab="Tests Conclusion", ylab="Substitutions")
boxplot(all_MLE$deltaSVM~all_MLE$diffTests, outline=F, notch=T, xlab="Tests Conclusion", ylab="deltaSVM")

################################################################################
### Boxplots
pdf(paste0(pathFigures, "/All_NarrowPeaks_boxplots_", method, "_add_4sub_treshold_", qval_treshold, ".pdf"), width = 10)
par(mfrow=c(1,2))
boxplot(all_MLE$Nmut~all_MLE$species, outline=F, xlab="", ylab="Substitutions", las=2)
boxplot(all_MLE$length~all_MLE$species, outline=F, xlab="", ylab="Length", las=2)

a <- boxplot(all_MLE$Prop.Sub~all_MLE$species, outline=F, xlab="", ylab="Substitution per base", las=2)
boxplot(all_MLE$Prop.Sub~all_MLE$TF, outline=F, xlab="TFs", ylab="Substitution per base", las=1)

boxplot(all_MLE$SVM~all_MLE$species, outline=F, xlab="", ylab="SVM", las=2)
boxplot(all_MLE$SVM~all_MLE$TF, outline=F, xlab="TFs", ylab="SVM", las=1)

boxplot(all_MLE$SumObs~all_MLE$species, outline=F, xlab="", ylab="Sum DeltaSVM Substitutions", las=2)
boxplot(all_MLE$SumObs~all_MLE$TF, outline=F, xlab="TFs", ylab="Sum DeltaSVM Substitutions", las=1)

boxplot(all_MLE$deltaSVM~all_MLE$species, outline=F, xlab="", ylab="DeltaSVM", las=2)
boxplot(all_MLE$deltaSVM~all_MLE$TF, outline=F, xlab="TFs", ylab="DeltaSVM", las=1)
dev.off()


Conclu.sp.TF <- table(all_MLE$TF, all_MLE$species, all_MLE$qval_Conclusion)
sp.TF <- table(all_MLE$TF, all_MLE$species)

Directional_count <- Conclu.sp.TF[,,1]+Conclu.sp.TF[,,2]
Stabilising_count <- Conclu.sp.TF[,,3]
Neutral_count <- Conclu.sp.TF[,,4]
prop_positive <- as.data.frame(Directional_count/sp.TF)

pdf(paste0(pathFigures, "/All_NarrowPeaks_Test_Conclusion_", method, "_qval0.1_4sub_treshold_", qval_treshold, ".pdf"), width = 10)
par(mfrow=c(1,2))
boxplot(prop_positive$Freq~prop_positive$Var2, las=3, ylab="Proportion of Directional signal", xlab="")
boxplot(prop_positive$Freq~prop_positive$Var1, las=3, ylab="Proportion of Directional signal", xlab="")

ggplot(prop_positive, aes(Var2, Var1, fill=Freq*100)) + 
  geom_tile(colour = "black") + 
  scale_fill_gradient(low="white", high="red", name = "% positive") +
  coord_equal(expand = 0) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        legend.text = element_text(size=15), legend.title = element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text=element_text(size=15), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Transcription factors", x = "Species")

prop_stab <- as.data.frame(Stabilising_count/sp.TF)
boxplot(prop_stab$Freq~prop_stab$Var2,las=3 , ylab="Proportion of Stabilising signal", xlab="")
boxplot(prop_stab$Freq~prop_stab$Var1,las=3, ylab="Proportion of Stabilising signal", xlab="")

prop_null <- as.data.frame(Neutral_count/sp.TF)
boxplot(prop_null$Freq~prop_null$Var2, las=3, ylab="Proportion of Neutral signal", xlab="")
boxplot(prop_null$Freq~prop_null$Var1, las=3, ylab="Proportion of Neutral signal", xlab="")

dev.off()

################################################################################
#### Bar plots
pdf(paste0(pathFigures, "/All_NarrowPeaks_models_comparisons_barplots_", method, "_qval0.1_4sub_treshold_", qval_treshold, ".pdf"), width = 11)
par(mfrow=c(1,2))

# Overlap Tests
conclu.test.Max <- prop.table(table(all_MLE$signif.FDR, all_MLE$qval_Conclusion))
barplot(conclu.test.Max, col=c("forestgreen", "#bae4b3", "firebrick"), las=1, ylab="Proportion", main="maxLL vs Permutation Test")

conclu.test <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$signif.FDR))
barplot(conclu.test, col=col, las=1, ylab="Proportion", main="Permutation vs maxLL Test")

# Conclusion * Species
Conclu.sp.permut <- prop.table(table(all_MLE$signif.FDR, all_MLE$exp), margin = 2)
ymax=max(Conclu.sp.permut[1,-12], na.rm=T)
barplot(Conclu.sp.permut[1:2,-12], col=col, las=2, ylab="Proportion", main="Permutation Test", ylim=c(0,ymax))
legend("topleft", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")

Conclu.sp <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$species), margin = 2)
barplot(Conclu.sp[1:2,-12], col=col, las=2, ylab="Proportion", main="MaxLL Test")

#alpha_treshold=0.1
#all_MLE$signif_FDR <- "Neutral"
#all_MLE$signif_FDR <- as.factor(ifelse(all_MLE$FDR.high <= alpha_treshold, "Directional", all_MLE$signif_FDR))
#all_MLE$signif_FDR <- as.factor(ifelse(all_MLE$FDR.low <= alpha_treshold, "Directional", all_MLE$signif_FDR))   
#all_MLE$signif_FDR <- factor(all_MLE$signif_FDR, levels=c("Directional", "Neutral"))

par(mfrow=c(2,2))
Conclu.sp <- prop.table(table(all_MLE$Conclusion, all_MLE$species), margin = 2)
barplot(Conclu.sp[1:2,-11], col=col, las=2, ylab="Proportion", main="MaxLL Test pval<0.01")

Conclu.sp <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$species), margin = 2)
barplot(Conclu.sp[1:2,-11], col=col, las=2, ylab="Proportion", main="MaxLL Test FDR<0.1")

Conclu.sp <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$species), margin = 2)
barplot(Conclu.sp[1:2,-11], col=col, las=2, ylab="Proportion", main="MaxLL Test qval<0.1")

#plot(all_MLE$qval_null_pos~all_MLE$FDR_null_pos)
#abline(a=0, b=1, col="red")

par(mfrow=c(2,2))
Conclu.sp.permut <- prop.table(table(all_MLE$signif, all_MLE$species), margin = 2)
ymax=max(Conclu.sp.permut[1,-12], na.rm=T)
barplot(Conclu.sp.permut[1:2,-12], col=col, las=2, ylab="Proportion", main="Permutation Test", ylim=c(0,ymax))

Conclu.sp.permut <- prop.table(table(all_MLE$signif.FDR, all_MLE$species), margin = 2)
ymax=max(Conclu.sp.permut[1,-12], na.rm=T)
barplot(Conclu.sp.permut[1,-12], col=col[1], las=2, ylab="Proportion", main="Permutation FDR<0.1", ylim=c(0,ymax))

Conclu.sp <- prop.table(table(all_MLE$signif.two, all_MLE$species), margin = 2)
barplot(Conclu.sp[1,-12], col="forestgreen", las=2, ylab="Proportion", main="Permut Two-tailed")

Conclu.sp <- prop.table(table(all_MLE$signif.two.FDR, all_MLE$species), margin = 2)
barplot(Conclu.sp[1,-12], col="forestgreen", las=2, ylab="Proportion", main="Permut FDR<0.1")

#plot(all_MLE$qval.high~all_MLE$FDR.high)
#abline(a=0, b=1, col="red")

# Conclusion * TFs
par(mfrow=c(2,2))
Conclu.TF.permut <- prop.table(table(all_MLE$signif, all_MLE$TF), margin = 2)
ymax=max(Conclu.TF.permut[1,], na.rm=T)+0.01
barplot(Conclu.TF.permut[1,], col="forestgreen", las=3, ylab="Proportion", main="Permutation Test", ylim=c(0,ymax))

Conclu.TF <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$TF), margin = 2)
barplot(Conclu.TF[1:2,], col=col, las=3, ylab="Proportion", main="MaxLL Test", ylim=c(0,ymax))
legend("topleft", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")

# Conclusion * Sample 
Conclu.samp.permut <- prop.table(table(all_MLE$signif, all_MLE$sample), margin = 2)
ymax=max(Conclu.samp.permut[1,], na.rm=T)
barplot(Conclu.samp.permut[1,], col="forestgreen", las=3, ylab="Proportion", main="Permutation Test", ylim=c(0,ymax))

Conclu.samp <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$sample), margin = 2)
barplot(Conclu.samp[1:2,], col=col, las=3, ylab="Proportion", main="MaxLL Test", ylim=c(0,ymax))
legend("topright", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")

# Conclusion * Substitutions
par(mfrow=c(2,2))
quantile_5 = paste0(seq(5, 100, 5), "%")
quantile_10 = paste0(seq(10, 100, 10), "%")
quantile_50 = paste0(seq(1, 100, 2), "%")
delta=expression(Delta)

all_MLE$Sub_class <- cut(all_MLE$Prop.Sub,  breaks=quantile(all_MLE$Prop.Sub, probs = seq(0, 1, 0.1)), include.lowest = T, dig.lab=1)

all_MLE2 = all_MLE

for (sp in species[species != "macaca"]){
  all_MLE = all_MLE2[which(all_MLE2$species==sp),]
  all_MLE$Sub_class <- cut(all_MLE$Prop.Sub,  breaks=quantile(all_MLE$Prop.Sub, probs = seq(0, 1, 0.02)), include.lowest = T, dig.lab=1)
  Conclu.Sub <- prop.table(table(all_MLE$Conclusion, all_MLE$Sub_class), margin = 2)
  Conclu.Sub.permut <- prop.table(table(all_MLE$signif, all_MLE$Sub_class), margin = 2)
  
  barplot(Conclu.Sub.permut[1:2,], col=col, las=1, ylab="Proportion", xlab="Substitution per base (quantile)",  main="Permutation Test", names=quantile_50)
  legend("topleft", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")
  
  barplot(Conclu.Sub[1:2,], col=col, las=1, ylab="Proportion", xlab="Substitution per base (quantile)",  main=sp, names=quantile_50)
  
}

all_MLE = all_MLE2
#barplot(Conclu.Sub[3,], col=col[3], las=1, ylab="Proportion", xlab="Substitution per base (quantile)",  main="MaxLL Test", names=quantile_10)
#barplot(Conclu.Sub[4,], col=col[4], las=1, ylab="Proportion", xlab="Substitution per base (quantile)", 
#        main="MaxLL Test", names=quantile_10, ylim=c(0,1))

all_MLE[which(all_MLE$Nmut >20),]$Nmut <- 20

#Conclusion * SumObsSVM
all_MLE$SumObs_class <- cut(all_MLE$SumObs,  breaks=quantile(all_MLE$SumObs, probs = seq(0, 1, 0.05)), include.lowest = T)

Conclu.SumObs.permut <- prop.table(table(all_MLE$signif, all_MLE$SumObs_class), margin = 2)
barplot(Conclu.SumObs.permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab=expression("Sum "*Delta*"SVM (quantile)"),
        main="Permutation Test", names=quantile_5)

Conclu.SumObs <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$SumObs_class), margin = 2)
barplot(Conclu.SumObs[1:3,], col=col, las=1, ylab="Proportion", xlab=expression("Sum "*Delta*"SVM (quantile)"), main="MaxLL Test")
legend("topleft", legend=c("Directional (+)", "Directional (-)", "Stabilising"), fill=col, bty="n")

# Conclusion * deltaSVM
all_MLE$deltaSVM_class <- cut(all_MLE$deltaSVM,  breaks=quantile(all_MLE$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)

Conclu.deltaSVM.permut <- prop.table(table(all_MLE$signif, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM.permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab="DeltaSVM (quantile)", main="Permutation Test", names=quantile_5)

Conclu.deltaSVM <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$deltaSVM_class), margin = 2)
barplot(Conclu.deltaSVM[1:2,], col=col, las=1, ylab="Proportion", xlab="DeltaSVM (quantile)", main="MaxLL Test")

# Conclusion * GC content
all_MLE$GC_class <- cut(all_MLE$GC_Content,  breaks=quantile(all_MLE$GC_Content, probs = seq(0, 1, 0.05)), include.lowest = T)
Conclu.GC_class_permut <- prop.table(table(all_MLE$signif, all_MLE$GC_class), margin = 2)
barplot(Conclu.GC_class_permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab="GC Content (quantile)", main="Permutation Test", names=quantile_5)
legend("topright", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")

Conclu.GC_class <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$GC_class), margin = 2)
barplot(Conclu.GC_class[1:2,], col=col, las=1, ylab="Proportion", xlab="GC Content (quantile)", main="MaxLL Test")

par(mfrow=c(1,2))
for (sub_type in c("Weak2Weak", "Weak2Strong", "Strong2Strong", "Strong2Weak")){
  print(sub_type)
  all_MLE[[paste0(sub_type, "_prop")]] <- all_MLE[[sub_type]]/all_MLE$Nmut
  sub_type = paste0(sub_type, "_prop")
  all_MLE$sub_class <- cut(all_MLE[[sub_type]],  breaks= seq(0,1,0.1), include.lowest = T)
  Conclu.sub_class <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$sub_class), margin = 2)
  barplot(Conclu.sub_class[1:2,], col=col, las=1, ylab="Proportion", xlab=paste("Quantile", sub_type), main="MaxLL Test")
  
  Conclu.sub_class_permut <- prop.table(table(all_MLE$signif, all_MLE$sub_class), margin = 2)
  barplot(Conclu.sub_class_permut[1,], col="forestgreen", las=1, ylab="Proportion", xlab=paste("Quantile", sub_type), main="Permutation Test")
}

# Epistatic effect
par(mfrow=c(1,2))
boxplot(all_MLE$Diffdelta~all_MLE$Sub_class, outline=F, xlab="Substitution per base (quantile)", ylab=expression("Additive - Total "*Delta*"SVM"), names=quantile_10)

all_MLE$Diffdelta_class <- cut(-all_MLE$Diffdelta,  breaks=c(min(all_MLE$Diffdelta), -1, -0.5, -0.25, 0, 0.25, 0.5, 1, max(all_MLE$Diffdelta)), include.lowest = T)
Conclu.epistatic <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$Diffdelta_class), margin = 2)
Conclu.epistatic.permut <- prop.table(table(all_MLE$signif, all_MLE$Diffdelta_class), margin = 2)
barplot(Conclu.epistatic[1:2,], col=col, las=2, ylab="Proportion", main="MaxLL Test")
legend("topright", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")
barplot(Conclu.epistatic.permut[1:2,], col=col, las=2, ylab="Proportion", main="Permut")
legend("topleft", legend=c("Directional (+)", "Directional (-)"), fill=col, bty="n")

plot(Conclu.epistatic.permut[1,]+Conclu.epistatic.permut[1,], type="b", col="red", ylim=c(0,0.2), las=1, 
     ylab="Proportion", xlab="Epistatic effect")
lines(Conclu.epistatic[1,]+Conclu.epistatic[1,], type="b")

dev.off()

## Lines
pdf(paste0(pathFigures, "/All_NarrowPeaks_comparisons_lines_", method, "_qval0.1_4sub_treshold_", qval_treshold, ".pdf"), width = 10)
par(mfrow=c(1,2))
pos = c(1, 5, 10, 15, 20)
#delta permut
plot(Conclu.SumObs.permut[1:1,], type="b", ylim=c(0,0.5), col="forestgreen", pch=19, las=1, ylab="Proportion", xlab=expression("Sum "*Delta*"SVM (quantile)"), axes=F)
axis(1, at=pos, labels=c("<10", "(-3,-2.2]", "(0.6,1.2]", "(4.5,5.7]", ">18"), las=1)
axis(2, las=2)

#delta
Directional <- Conclu.SumObs[1,]+Conclu.SumObs[2,]
plot(Directional, type="b", ylim=c(0,0.5), col=col[1], pch=19, las=1, ylab="Proportion",
     xlab=expression("Sum "*Delta*"SVM (quantile)"), axes=F, main="")
lines(Conclu.SumObs[3,], type="b", col=col[3], pch=19)
axis(1, at=pos, labels=c("<10", "(-3,-2.2]", "(0.6,1.2]", "(4.5,5.7]", ">18"), las=1)
axis(2, las=2)
legend("topleft", legend=c("Directional", "Stabilising"), fill=c(col[1],col[3]), bty="n")

#Sub permut
plot(Conclu.Sub.permut[1:1,], type="b", ylim=c(0,0.13), col=col[1], pch=19, las=1, ylab="Proportion", xlab="Substitution per base (quantile)", axes=F)
axis(1, at=pos, labels=c("<5e-3", "(8e-3,9e-3]", "(1.5e-2,1.6e-2]", "(3e-2,4e-2]", ">9e-2"), las=1)
axis(2, las=2)

#Sub
Directional <- Conclu.Sub[1,]+Conclu.Sub[2,]
plot(Directional, type="b", ylim=c(0,0.13), col=col[1], pch=19, las=1, ylab="Proportion", xlab="Substitution per base (quantile)", axes=F)
lines(Conclu.Sub[3,], type="b", col=col[3], pch=19)
axis(1, at=pos, labels=c("<5e-3", "(8e-3,9e-3]", "(1.5e-2,1.6e-2]", "(3e-2,4e-2]", ">9e-2"), las=1)
axis(2, las=2)
legend("topleft", legend=c("Directional", "Stabilising"), col=c(col[1],col[3]), bty="n", lty=1, lwd=3)

# SVM
all_MLE$SVM_class <- cut(all_MLE$SVM,  breaks=quantile(all_MLE$SVM, probs = seq(0, 1, 0.05)), include.lowest = T)
Conclu.SVM_class_permut <- prop.table(table(all_MLE$signif.FDR, all_MLE$SVM_class), margin = 2)
Conclu.SVM_class <- prop.table(table(all_MLE$qval_Conclusion, all_MLE$SVM_class), margin = 2)
# Permut
plot(Conclu.SVM_class_permut[1:1,], type="b", ylim=c(0,0.1), col="forestgreen", pch=19, las=1, ylab="Proportion", xlab="SVM", axes=F)
axis(1, at=pos, labels=c("<265", "(-148, 130]", "(-85, -76]", "(-47, -41]", ">-4"), las=1)
axis(2, las=2)

# MaxLL
Directional <- Conclu.SVM_class[1,]+Conclu.SVM_class[2,]
plot(Directional, type="b", ylim=c(0,0.04), col=col[1], pch=19, las=1, ylab="Proportion", xlab="SVM", axes=F)
lines(Conclu.SVM_class[3,], type="b", col=col[3], pch=19)
axis(1, at=pos, labels=c("<265", "(-148, 130]", "(-85, -76]", "(-47, -41]", ">-4"), las=1)
axis(2, las=2)
legend("topleft", legend=c("Directional", "Stabilising"), col=c(col[1],col[3]), bty="n", lty=1, lwd=3)

# GC permut
plot(Conclu.GC_class_permut[1:1,], type="b", ylim=c(0,0.08), col="forestgreen", pch=19, las=1, ylab="Proportion", xlab="GC Content (%)", axes=F)
axis(1, at=pos, labels=c("<35", "(40, 42]", "(45, 47]", "(50, 52]", ">60"), las=1)
axis(2, las=2)
# GC
Directional <- Conclu.GC_class[1,]+Conclu.GC_class[2,]
plot(Directional, type="b", ylim=c(0,0.08), col=col[1], pch=19, las=1, ylab="Proportion", xlab="GC Content (%)", axes=F)
lines(Conclu.GC_class[3,], type="b", col=col[3], pch=19)
axis(1, at=pos, labels=c("<35", "(40, 42]", "(45, 47]", "(50, 52]", ">60"), las=1)
axis(2, las=2)
legend("topleft", legend=c("Directional", "Stabilising"), col=c(col[1],col[3]), bty="n", lty=1, lwd=3)

dev.off()
 
################################################################################
###### Analyse FDR and qvalues
pdf(paste0(pathFigures, "/qqplot_stab_pos_all_sp_4sub_treshold_", qval_treshold, ".pdf"), width=10)
par(mfrow=c(2,4))
for (sp in species){
  print(sp)
  for (TF in TFs){
    bigdata = all_MLE[which(all_MLE$species==sp & all_MLE$TF==TF),]
    #data$Sub_class <- cut(bigdata$Prop.Sub,  breaks=quantile(bigdata$Prop.Sub, probs = seq(0, 1, 0.1)), include.lowest = T, dig.lab=1)
    
    for (sub in c(2, 5, 7)){ #levels(as.factor(bigdata$Nmut)
      data = bigdata[which(bigdata$Nmut>=as.numeric(sub)),]
      if (nrow(data)>100){
        # Permut
        #qval= qvalue(data$pval.high)$qvalues
        #FDR = p.adjust(data$pval.high, method="fdr")
        #hist(data$pval.high, breaks=50,main="", xlab="pval. Permut")
        #hist(qval, breaks=20, main=paste(sp, TF, "Permut one-tailed"), xlab="qval. Permut", xlim=c(0,1))
        #hist(FDR, breaks=50, main="", xlab="FDR Permut", xlim=c(0,1))
        
        # Permut two-tailed
        #qval= qvalue(data$pval.two.tailed)$qvalues
        #FDR = p.adjust(data$pval.two.tailed, method="fdr")
        #hist(data$pval.two.tailed, breaks=50, main="", xlab="pval. Permut")
        #hist(qval, breaks=50, main=paste(sp, TF, "Permut two-tailed"), xlab="qval. Permut", xlim=c(0,1))
        #hist(FDR, breaks=50, main="", xlab="FDR Permut", xlim=c(0,1))
        
        # MLL Purif-pos
        #qval= qvalue(data$p_value_purif_pos)$qvalues
        #FDR = p.adjust(data$p_value_purif_pos, method="fdr")
        #hist(data$p_value_purif_pos, breaks=100, main=paste0(sp, TF, ", Nsub>=",sub,", N=", nrow(data)), xlab="pval. MLL Stab-Pos")
        qqunif(data$p_value_purif_pos, ci=T, main=paste0(sp, TF, ", Nsub>=",sub,", N=", nrow(data)))
        #signif = length(which(qval<0.2))
        #prop=signif(signif/length(qval), digits=3)*100
        #hist(qval, breaks=50, main=paste0("Nsignif=",signif," (",prop, "%)"), xlab="qval. MLL Stab-Pos", xlim=c(0,1))
        #abline(v=0.2, col="red")
        #hist(FDR, breaks=50, main="", xlab="FDR MLL Purif-Pos", xlim=c(0,1))
        
        # MLL Null-pos
        #qval= qvalue(data$p_value_null_pos)$qvalues
        #FDR = p.adjust(data$p_value_null_pos, method="fdr")
        #hist(data$p_value_null_pos, breaks=100, main="", xlab="pval. MLL Null-Pos")
        #hist(qval, breaks=50, main=paste(sp, TF, "Test MLL"), xlab="qval. MLL Null-Pos", xlim=c(0,1))
        #hist(FDR, breaks=50, main="", xlab="FDR MLL Null-Pos", xlim=c(0,1))
        
        # MLL Purif-Null
        #qval= qvalue(data$p_value_null_purif)$qvalues
        #FDR = p.adjust(data$p_value_null_pos, method="fdr")
        #hist(data$p_value_null_purif, breaks=100, main="", xlab="pval.ƒ MLL Null-Purif")
        #hist(qval, breaks=50, main=paste(sp, TF, "Test MLL"), xlab="qval. MLL Null-Purif", xlim=c(0,1))
        #hist(FDR, breaks=50, main="", xlab="FDR MLL Null-Pos", xlim=c(0,1))
        
      }
  }
}
}

dev.off()

################################################################################
###### Analyse non-similar conclusion
#par(mfrow=c(2,2))
#for (sp in species[species != "macaca"]){
#  subdata = all_MLE[which(all_MLE$species == sp & all_MLE$Nmut <21),]
#  cor <- cor.test(subdata$deltaSVM, subdata$Prop.Sub)$estimate
#  boxplot(subdata$deltaSVM~subdata$Nmut, main=sp, outline=F, xlab="Substitutions", ylab="DeltaSVM")
#  mtext(paste("R=", round(cor, 2)), side=3)
#}

#plot(Conclu.Sub.permut[1,]+Conclu.Sub.permut[2,], type="b", col="red", ylim=c(0,0.7), las=1, 
#     ylab="Proportion", xlab="Substitution per base (quantile)")
#lines(Conclu.Sub[1,]+Conclu.Sub[2,], type="b")

#all_MLE$similar <- ifelse(all_MLE$Conclusion==all_MLE$signif.FDR, 1, 0)

#par(mfrow=c(2,2))
# Proportion of similar per Sub_class
#subdataset <- all_MLE[which(all_MLE$diffTests != "Null/Null"),]
#subdataset <- subdataset[which(subdataset$diffTests != "Dir/Dir"),]
#subdataset <- subdataset[which(subdataset$Conclusion != "Stabilizing"),]
#similar_Sub <- prop.table(table(subdataset$similar, subdataset$Sub_class), margin = 2)
#similar_delta <- prop.table(table(subdataset$similar, subdataset$deltaSVM_class), margin = 2)
#similar_svm <- prop.table(table(subdataset$similar, subdataset$SVM_class), margin = 2)
#similar_epistasis <- prop.table(table(subdataset$similar, subdataset$Diffdelta_class), margin = 2)
#plot(similar_svm[2,], type="b", ylim=c(0,1), las=1, ylab="Proportion similar", xlab="SVM (quantile)", main="Both non-Neutral")
#plot(similar_Sub[2,], type="b", ylim=c(0,1), las=1, ylab="Proportion similar", xlab="deltaSVM (quantile)", main="Both non-Neutral")
#plot(similar_delta[2,], type="b", ylim=c(0,1), las=1, ylab="Proportion similar", xlab="Substitution (quantile)", main="Both non-Neutral")
#plot(similar_epistasis[2,], type="b", ylim=c(0,1), las=1, ylab="Proportion similar", xlab="Epistasis (quantile)", main="Both non-Neutral")

#subdataset <- subdataset[which(subdataset$deltaSVM_class %in% levels(subdataset$deltaSVM_class)[19:20]),]
#pdf(paste0(pathFigures, "/Stats_different_conclusion_allPeaks_4sub_treshold_", qval_treshold, ".pdf"))
#par(mfrow=c(2,2))
#sub_permut <- prop.table(table(subdataset$signif.FDR, subdataset$Sub_class), margin = 2)
#sub_MLL <- prop.table(table(subdataset$Conclusion, subdataset$Sub_class), margin = 2)
#plot(sub_permut[1,]+sub_permut[2,], type="b", col="red", ylim=c(0,1), las=1, ylab="Proportion Dir", xlab="Substitution (quantile)")
#lines(sub_MLL[1,]+sub_MLL[2,], type="b")

#delta_MLL <- prop.table(table(subdataset$Conclusion, subdataset$deltaSVM_class), margin = 2)
#plot(delta_permut[1,]+delta_permut[2,], type="b", col="red", ylim=c(0,1), las=1, ylab="Proportion Dir", xlab="deltaSVM (quantile)")
#delta_permut <- prop.table(table(subdataset$signif.FDR, subdataset$deltaSVM_class), margin = 2)
#lines(delta_MLL[1,]+delta_MLL[2,], type="b")
#boxplot(subdataset$deltaSVM~subdataset$Sub_class, outline=F, xlab="Substitution per base (quantile)")

#sub_permut <- prop.table(table(subdataset$signif.FDR, subdataset$Diffdelta_class), margin = 2)
#sub_MLL <- prop.table(table(subdataset$Conclusion, subdataset$Diffdelta_class), margin = 2)
#     ylab="Proportion Dir", xlab="Epistasis (quantile)")
#plot(sub_permut[1,]+sub_permut[2,], type="b", col="red", ylim=c(0,1), las=1, main="Synergetic to Compensatory",
#lines(sub_MLL[1,]+sub_MLL[2,], type="b")

#sub_permut <- prop.table(table(subdataset$signif.FDR, subdataset$GC_class), margin = 2)
#sub_MLL <- prop.table(table(subdataset$Conclusion, subdataset$GC_class), margin = 2)
#plot(sub_permut[1,]+sub_permut[2,], type="b", col="red", ylim=c(0,1), las=1, main="Low to High GC%", ylab="Proportion Dir", xlab="GC (quantile)")
#lines(sub_MLL[1,]+sub_MLL[2,], type="b")
#dev.off()
