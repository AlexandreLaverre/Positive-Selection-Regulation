
path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
sp <- "human"
TF <- "Wilson/CEBPA"

file = paste0(path, "positive_selection/NarrowPeaks/", sp, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
fileMaxLL = paste0(path, "positive_selection/NarrowPeaks/", sp, "/", TF, "/Tests/MLE_summary_exact.csv")

rand <- read.table(file, h=T)
rownames(rand) <- rand$ID
maxLL <- read.csv(fileMaxLL, h=T, row.names = 1)
maxLL$binary <- ifelse(grepl("Directional", maxLL$Conclusion), 1, 0) 

rand$MaxLL <- maxLL[rownames(rand),]$binary 
rand$classSVM <- cut(rand$SVM, breaks=quantile(rand$SVM, probs = seq(0, 1, 0.1)), include.lowest = T)
rand$signif.bin <- ifelse(rand$pval.high <= 0.01, 1, 0)
rand$positiveDelta <- ifelse(rand$deltaSVM > 0, 1, 0)
rand$FDR <- p.adjust(rand$pval.high, method="fdr")
rand$signif.FDR.bin <- ifelse(rand$FDR <= 0.01 , 1, 0)
rand <- rand %>% mutate(length = sapply(strsplit(ID, ":"), function(x) as.numeric(x[5]) - as.numeric(x[4])))
rand$classlength <- cut(rand$length, breaks=quantile(rand$length, probs = seq(0, 1, 0.1)), include.lowest = T)

a <- table(rand$classlength,rand$signif.bin)
barplot(a[,2]*100/(a[,2]+a[,1]), las=1, xlab="Length quantile", ylab="% signif")

a <- table(rand$classSVM,rand$signif.bin)
barplot(a[,2]*100/(a[,2]+a[,1]), las=1, xlab="SVM quantile", ylab="% signif")

a <- table(rand$classSVM,rand$signif.FDR.bin)
barplot(a[,2]*100/(a[,2]+a[,1]), las=1, xlab="SVM quantile", ylab="% signif FDR")

all_prop_signif <- c()
all_prop_positif <- c()
for (q in seq(0,0.9,0.1)){
  sub <- rand[which(rand$SVM >= quantile(rand$SVM, q)),]
  prop_signif <- sum(sub$signif.bin)/nrow(sub)
  prop_positif <-  sum(sub$positiveDelta)/nrow(sub)
  
  all_prop_signif <- c(all_prop_signif, prop_signif)
  all_prop_positif <- c(all_prop_positif, prop_positif)
}

par(mfrow=c(1,2))
#barplot(all_prop_positif, las=1, names=paste0(100-seq(0, 90, 10), "%"), cex.names=0.75,  ylim=c(0, 0.7),
#        xlab="Fraction of dataset (highest SVM)", ylab="Proportion positive deltaSVM")
barplot(all_prop_signif, las=1, names=paste0(100-seq(0, 90, 10), "%"), cex.names=0.75, ylim=c(0, 0.1),
        xlab="Fraction of dataset (highest SVM)", ylab="Proportion significant peaks", main="Permutation Test")


### MaxLL Test
all_prop_signif <- c()
all_prop_positif <- c()
for (q in seq(0,0.9,0.1)){
  sub <- rand[which(rand$SVM >= quantile(rand$SVM, q)),]
  prop_signif <- sum(sub$MaxLL)/nrow(sub)
  prop_positif <-  sum(sub$positiveDelta)/nrow(sub)
  
  all_prop_signif <- c(all_prop_signif, prop_signif)
  all_prop_positif <- c(all_prop_positif, prop_positif)
}


barplot(all_prop_signif, las=1, names=paste0(100-seq(0, 90, 10), "%"), cex.names=0.75, ylim=c(0, 0.1),
        xlab="Fraction of dataset (highest SVM)", ylab="Proportion significant peaks", main="MaxLL Test")

#dev.off()