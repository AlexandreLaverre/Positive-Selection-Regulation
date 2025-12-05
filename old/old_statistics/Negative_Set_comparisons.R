library(corrplot)
library(RColorBrewer)
library(stringr)

par(xpd = TRUE)
path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Neutral model")
col_permut=c("forestgreen", "deepskyblue3", "firebrick")
col_LRT=c("firebrick", "forestgreen", "white", "deepskyblue3")
names(col_permut) = c("Pos Increase", "Pos Decrease", "Neutral")
maxSub = 150

################################################################################
####### Data
obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
deltas <- read.table(paste0(path, "/positive_selection/BroadPeaks/human/delta_negative_set/CEBPA/deltas/focal_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
row.names(deltas) <- deltas$seq_name
  

#### Permutation Test
permut <- read.table(paste0(path, "positive_selection/BroadPeaks/human/delta_negative_set/CEBPA/PosSelTest_deltaSVM_10000permutations.txt"), h=T, row.names = 1)
permut$signif <- as.factor(ifelse(permut$pval.high <= 0.01, "Pos Increase", ifelse(permut$pval.high >= 0.99, "Pos Decrease", "Neutral")))
permut$signif.bin <- ifelse(permut$pval.high <= 0.01, 1, ifelse(permut$pval.high >= 0.99, -1, 0))
permut$signif <- factor(permut$signif, levels=c("Pos Increase", "Pos Decrease", "Neutral"))

#### MaxLL
negative <- read.csv(paste0(path, "MaxLikelihoodApproach/human/CEBPA/simulations/MLE_summary_negative_set_50bins.csv"), h=T, row.names = 1)
common <- intersect(rownames(permut), rownames(negative))
negative <- negative[common,]
permut <- permut[common, ]

negative$Conclusion <- as.factor(negative$Conclusion)
negative$Conclusion <- factor(negative$Conclusion, levels=c("Positive model", "Positive model Decreasing", "Positive model Increasing", "Stabilizing model", "Neutral model"))
negative$ID <- rownames(negative)
negative$length <- as.numeric(str_split_i(negative$ID, ":", 3))-as.numeric(str_split_i(negative$ID, ":", 2))
negative$deltaSVM <- deltas[row.names(negative),]$deltaSVM
negative$Diffdelta <-negative$deltaSVM-negative$SumObs

negative <- cbind(negative, permut[,c("SVM", "med.expected.deltaSVM", "mean.expected.deltaSVM", "pval.high", "signif", "signif.bin")])
negative$classSVM <- cut(negative$SVM, breaks=quantile(negative$SVM, probs = seq(0, 1, 0.1)), include.lowest = T)
negative$classdeltaSVM <- cut(negative$deltaSVM, breaks=quantile(negative$deltaSVM, probs = seq(0, 1, 0.1)), include.lowest = T)
negative$classlength <- cut(negative$length, breaks=quantile(negative$length, probs = seq(0, 1, 0.1)), include.lowest = T)
negative$classMut <- cut(negative$Nmut, breaks=c(0,2,3,4,5,8,10,15, max(negative$Nmut)), include.lowest = T)

negative$Model.bin <- ifelse(negative$Conclusion == "Positive model", 1, 0)
negative$positiveDelta <- ifelse(negative$deltaSVM > 0, 1, 0)
#rownames(negative) <- NULL

negative$direction <- as.factor(ifelse(negative$SumObs>negative$med.expected.deltaSVM, "Increasing", "Decreasing"))
negative$Combine <- paste(negative$Conclusion, negative$direction)
negative[which(negative$Conclusion=="Positive model"),]$Conclusion <- negative[which(negative$Conclusion=="Positive model"),]$Combine
negative$Conclusion <- factor(negative$Conclusion, levels=c("Positive model Decreasing", "Positive model Increasing", "Stabilizing model", "Neutral model"))

print(table(negative$Conclusion))
print(table(negative$signif))

pdf(paste0(path, "figures/Negative_Set_Comparisons_human_CEBPA.pdf"), width = 10)
################################################################################
#### Likelihood Ratio Test summary
par(mfrow=c(2,2))
barplot(table(negative$Conclusion)/nrow(negative), col=col_LRT, las=1, main="Negative set",
        xlab="LRT Conclusion", ylab="Proportion", cex.lab=1.2, cex.names=1.2)
boxplot(negative$MeanObs~negative$Conclusion, notch=T, outline=F, col=col_LRT,
        xlab="LRT conclusion", ylab="Mean Delta Obs",cex.lab=1.2, cex.names=1.2)
boxplot(negative$VarObs~negative$Conclusion, notch=F, outline=F, col=col_LRT,
        xlab="LRT conclusion", ylab="Var Delta Obs", cex.lab=1.2, cex.names=1.2)
boxplot(negative$Nmut~negative$Conclusion, notch=F, outline=F, col=col_LRT,
        xlab="LRT conclusion", ylab="Nb subsitution",cex.lab=1.2, cex.names=1.2)

##### Permutations test summary
par(mfrow=c(2,2))
barplot(table(negative$signif)/nrow(negative), col=col_permut, las=1, main="Negative set",
        xlab="Permut Conclusion", ylab="Proportion", cex.lab=1.2, cex.names=1.2)
boxplot(negative$MeanObs~negative$signif, notch=F, outline=F, col=col_permut,
        xlab="Permut conclusion", ylab="Mean Delta Obs", cex.lab=1.2, cex.names=1.2)
boxplot(negative$VarObs~negative$signif, notch=F, outline=F, col=col_permut,
        xlab="Permut conclusion", ylab="Var Delta Obs", cex.lab=1.2, cex.names=1.2)
boxplot(negative$Nmut~negative$signif, notch=F, outline=F, col=col_permut,
        xlab="Permut conclusion", ylab="Nb subsitution",cex.lab=1.2, cex.names=1.2)

################################################################################
#### Impact mutations
par(mfrow=c(1,2))
a <- table(negative$Nmut, negative$signif.bin)
barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="Nb substitutions", 
        ylab="Proportion signif", main="Permut Test")

a <- table(negative$Nmut, negative$Model.bin)
barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="Nb substitutions", 
        ylab="Proportion signif", main="LL Ratio Test")


######### Permutations ######### 
#### Impact SVM
par(mfrow=c(2,2))
par(xpd = TRUE)
col_permut=c("firebrick", "deepskyblue3", "forestgreen")
a <- table(negative$classSVM, negative$signif.bin)/nrow(negative)*10
barplot(t(a), col=col_permut,
        las=1, xlab="SVM quantile", ylab="Proportion", main="Permutations Test")
legend("top", legend=c("Positive (-)", "Neutral", "Positive (+)"), fill=col_permut, bty="n", ncol=3, inset = c(0, -0.15), cex=1)

a <- table(negative$classdeltaSVM, negative$signif.bin)/nrow(negative)*10
barplot(t(a), col=col_permut,
        las=1, xlab="deltaSVM quantile", ylab="Proportion", main="Permutations Test")

####  Length
a <- table(negative$classlength, negative$signif.bin)/as.numeric(table(negative$classlength))
barplot(t(a), col=col_permut,
        las=1, xlab="Length quantile", ylab="Proportion", main="Permutations Test")
legend("top", legend=c("Positive (-)", "Neutral", "Positive (+)"), fill=col_permut, bty="n", ncol=3, inset = c(0, -0.15), cex=1)

# Mutations
a <- table(negative$classMut, negative$signif.bin)/as.numeric(table(negative$classMut))
barplot(t(a), col=col_permut,
        las=1, xlab="Nb substitutions", ylab="Proportion", main="Permutations Test")

######### LRT ######### 
par(mfrow=c(2,2))
par(xpd = TRUE)

a <- table(negative$classSVM, negative$Conclusion)/nrow(negative)*10
barplot(t(a), col=col_LRT, las=1, xlab="SVM quantile", 
        ylab="Proportion", main="LL Ratio Test")
legend("top", legend=c("Positive (-)", "Positive (+)", "Stabilising", "Neutral"), fill=col_LRT, bty="n", ncol=4, inset = c(0, -0.15), cex=1)

a <- table(negative$classdeltaSVM, negative$Conclusion)/nrow(negative)*10
barplot(t(a), col=col_LRT, las=1, xlab="deltaSVM quantile", ylab="Proportion", main="LL Ratio Test")

a <- table(negative$classlength, negative$Conclusion)/as.numeric(table(negative$classlength))
barplot(t(a), col=col_LRT,las=1, xlab="Length quantile", ylab="Proportion", main="LL Ratio Test")
legend("top", legend=c("Positive (-)", "Positive (+)", "Stabilising", "Neutral"), fill=col_LRT, bty="n", ncol=4, inset = c(0, -0.15), cex=1)

a <- table(negative$classMut, negative$Conclusion)/as.numeric(table(negative$classMut))
barplot(t(a), col=col_LRT,las=1, xlab="Nb substitutions", ylab="Proportion", main="LL Ratio Test")


################################################################################
par(mfrow=c(1,2))
barplot(table(negative$direction, negative$Conclusion)/nrow(negative), 
        col=c("firebrick", "forestgreen"), ylab="Proportion", main="LL Ratio Test")
legend("topleft", legend=c("Decreasing", "Increasing"), fill=c("firebrick", "forestgreen"), bty="n")

barplot(table(negative$direction, negative$signif)/nrow(negative), 
        col=c("firebrick", "forestgreen"), ylab="Proportion", main="Permut Test")
dev.off()