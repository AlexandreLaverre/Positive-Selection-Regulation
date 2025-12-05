library(stringi)

# check distrib delta and correlation
sp = "human"
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/"
positive <- read.table(paste0(path, sp, "/PosSel/PosSelTest_deltaSVM_10000permutations.txt"), h=T)
positive$FDR <- p.adjust(positive$pval.high, method="fdr")

par(mfrow=c(2,2))
par(mar=c(4.5,4.5,1,1))
hist(positive$SVM, breaks=100, xlab="SVM score", main="human CEBPA positive")
hist(positive$deltaSVM, breaks=100, xlab="delta SVM", xlim=c(-20,20), main="")
#plot(positive$deltaSVM~positive$SVM, cex=0.5, xlab="SVM score", ylab="deltaSVM")
#R2 = cor.test(positive$deltaSVM, positive$SVM)$estimate
#mtext(paste0("R2=", signif(R2,2)), side=1, line=-2)

positive$signif <- ifelse(positive$FDR<0.1, "red", "black")
plot(positive$deltaSVM~positive$pval.high, cex=0.5, xlab="pval high", ylab="deltaSVM", col=positive$signif)
legend("topright", "FDR<0.1", col="red", pch=1, bty="n")
R2 = cor.test(positive$deltaSVM, positive$pval.high)$estimate
signif = nrow(positive[which(positive$FDR<0.1),])
mtext(paste0("R2=", signif(R2,2), "; ", "Nsignif=", signif), side=1, line=-2)


plot(positive$med.deltaSVM.simul~positive$SVM, ylim=c(-50, 1), xlab="SVM score", ylab="deltaSVM simul")

positive$signif <- ifelse(positive$pval.high<0.05, "green", "black")
plot(log(positive$deltaSVM)~positive$NbSub, xlim=c(0, 50), cex=0.4, col=positive$signif)
cor.test(positive$deltaSVM, positive$NbSub)

plot(positive$pval.high~positive$NbSub, xlim=c(0, 50), cex=0.4)
cor.test(positive$pval.high, positive$NbSub)

positive$class_SVM <- cut(positive$SVM, breaks=c(-1000, -400, -200, -100, 100), labels=c( "very_low", "low", "medium", "high"), include.lowest = T)
positive$class_SVM <- cut2(positive$SVM, g=4, labels=c( "very_low", "low", "medium", "high"), include.lowest = T)

# Length
positive$length <- as.numeric(str_split_i(positive$ID, ":", 3))-as.numeric(str_split_i(positive$ID, ":", 2))
positive$Sub_ratio <- positive$NbSub/positive$length
plot(positive$SVM~positive$length, cex=0.4, xlab="Seq length", ylab="SVM score")
plot(positive$deltaSVM~positive$length, cex=0.4, xlab="Seq length", ylab="delta SVM")
cor.test(positive$deltaSVM, positive$length)
cor.test(positive$pval.high, positive$length)
plot(positive$NbSub~positive$length, cex=0.4, xlab="Seq length", ylab="Nb substitutions", ylim=c(0,50))
cor.test(positive$NbSub, positive$length)


negative <- read.table(paste0(path, sp, "/delta_negative_set/CEBPA/PosSelTest_deltaSVM_1000permutations.txt"), h=T)
negative$FDR <- p.adjust(negative$pval.high, method="fdr")

par(mfrow=c(2,2))
boxplot(positive$SVM, negative$SVM, notch=T, outline=F, names=c("Positive set", "Negative set"), ylab="SVM")
boxplot(positive$deltaSVM, negative$deltaSVM, notch=T, outline=F, names=c("Positive set", "Negative set"), ylab="deltaSVM")
#hist(negative$SVM, breaks=100, xlab="SVM score", main="human CEBPA negative")
#hist(negative$deltaSVM, breaks=100, xlab="delta SVM", xlim=c(-20,20), main="")
#plot(negative$deltaSVM~negative$SVM, cex=0.5, xlab="SVM score", ylab="deltaSVM", main="Human CEBPA Negative")
#R2 = cor.test(negative$deltaSVM, negative$SVM)$estimate
#mtext(paste0("R2=", signif(R2,2)), side=1, line=-2)
hist(negative$deltaSVM, breaks=100, xlab="delta SVM", xlim=c(-20,20), main="Human CEBPA Negative")
mtext(paste0("median=", signif(median(negative$deltaSVM),2), "\n", "mean=", signif(mean(negative$deltaSVM),2)), at=-12, line=-3, cex=0.8)

negative$signif <- ifelse(negative$FDR<0.1, "red", "black")
plot(negative$deltaSVM~negative$pval.high, cex=0.5, xlab="pval high", ylab="deltaSVM", col=negative$signif)
legend("topright", "FDR<0.1", col="red", pch=1, bty="n")
R2 = cor.test(negative$deltaSVM, negative$pval.high)$estimate
signif = nrow(negative[which(negative$FDR<0.1),])
mtext(paste0("R2=", signif(R2,2), "; ", "Nsignif=", signif), side=1, line=-2)


plot(negative$deltaSVM~negative$pval.high, cex=0.4)
plot(negative$pval.high~negative$NbSub, xlim=c(0, 50), cex=0.4)
cor.test(negative$pval.high, negative$NbSub)

plot(negative$deltaSVM~negative$NbSub, xlim=c(0, 50), cex=0.4)
cor.test(negative$pval.high, negative$NbSub)
# test score sequence via gkmpredict
## Positive 
positive_score <- read.table(paste0(path, sp, "/CEBPA/sequences/positive_SVMscore.txt"), h=F)
colnames(positive_score) <- c("ID", "score")
rownames(positive_score) <- positive_score$ID
positive_score <- positive_score[positive$ID,]

plot(positive$SVM~positive_score$score, cex=0.4, ylab="sum kmer SVM", xlab="predict SVM", main = "Positive")
cor.test(positive$SVM, positive_score$score)

## Ancestral  
ancestral_score <- read.table(paste0(path, sp, "/CEBPA/sequences/positive_ancestral_SVMscore.txt"), h=F)
colnames(ancestral_score) <- c("ID", "score")
rownames(ancestral_score) <- ancestral_score$ID
ancestral_score <- ancestral_score[positive$ID,]
ancestral_score$SVM <- positive$SVM + positive$deltaSVM
  
plot(ancestral_score$SVM~ancestral_score$score, cex=0.4, ylab="sum kmer SVM", xlab="predict SVM", main = "Ancestral")
cor.test(ancestral_score$SVM, ancestral_score$score)

## Negative  
negative_score <- read.table(paste0(path, sp, "/delta_negative_set/CEBPA/sequences/negative_SVMscore.txt"), h=F)
colnames(negative_score) <- c("ID", "score")
rownames(negative_score) <- negative_score$ID
negative_score <- negative_score[negative$ID,]

plot(negative$SVM~negative_score$score, cex=0.4, ylab="sum kmer SVM", xlab="predict SVM", main = "Negative")
cor.test(negative$SVM,negative_score$score)

boxplot(positive$SVM, ancestral_score$SVM, negative$SVM, notch=T, outline=F, ylab="sum kmer SVM",
        col=c("forestgreen", "lightblue3", "firebrick"), names=c("Positive", "Ancestral", "Negative"))

boxplot(positive_score$score, ancestral_score$score, negative_score$score, notch=T, outline=F, ylab="predict SVM",
        col=c("forestgreen", "lightblue3", "firebrick"), names=c("Positive", "Ancestral", "Negative"))

boxplot(positive$deltaSVM,negative$deltaSVM, notch=T, outline=F, ylab="deltaSVM",
        names=c("Positive", "Negative"), col=c("forestgreen", "firebrick"))

## Comparison deltaSVM versus delta predictded score (gkmpredict)
positive$delta_score <- positive_score$score-ancestral_score$score

plot(positive$delta_score~positive$deltaSVM, cex=0.5, xlab="deltaSVM", ylab="delta score")
cor.test(positive$delta_score,positive$deltaSVM)

positive$quali_deltaSVM <- ifelse(positive$deltaSVM > 0, "deltaSVM_positive", "deltaSVM_negative")
positive$quali_deltascore <- ifelse(positive$delta_score > 0, "delta_score_positive", "delta_score_negative")
table(positive$quali_deltaSVM, positive$quali_deltascore)


## Correlation GC ratio
pos_seq <- readDNAStringSet(paste0(path, sp, "/CEBPA/sequences/filtered_focal_sequences.fa"))

# Calculate the GC ratio for each sequence
gc_ratios <- letterFrequency(pos_seq, letters = c("G", "C")) / width(pos_seq)
gc_ratios <- data.frame(SequenceID = names(pos_seq), GCRatio = gc_ratios)
gc_ratios$GC_ratio = gc_ratios$GCRatio.G+gc_ratios$GCRatio.C
rownames(gc_ratios) <- gc_ratios$SequenceID
positive$GC_ratio <- gc_ratios[positive$ID,]$GC_ratio

plot(positive$SVM~positive$GC_ratio, cex=0.5)
cor.test(positive$SVM, positive$GC_ratio) # non signif
hist(positive$GC_ratio)
plot(positive$deltaSVM~positive$GC_ratio, cex=0.5) 
cor.test(positive$deltaSVM, positive$GC_ratio) # non signif

plot(positive$GC_ratio~positive$Sub_ratio, cex=0.5)
cor.test(positive$GC_ratio, positive$Sub_ratio) # non signif


## Substitutions
positive_substi <- read.table(paste0(path, sp, "/CEBPA/sequences/substitutions.txt"), h=T)
negative_substi <- read.table(paste0(path, sp, "/delta_negative_set/CEBPA/sequences/substitutions.txt"), h=T)
downstream_substi <- read.table(paste0(path, sp, "/CEBPA_flanking/downstream/sequences/substitutions.txt"), h=T)
upstream_substi <- read.table(paste0(path, sp, "/CEBPA_flanking/upstream/sequences/substitutions.txt"), h=T)

positive_substi[,4:7] <- positive_substi[,4:7]/positive_substi$Length
negative_substi[,4:7] <- negative_substi[,4:7]/negative_substi$Length
downstream_substi[,4:7] <- downstream_substi[,4:7]/downstream_substi$Length
upstream_substi[,4:7] <- upstream_substi[,4:7]/upstream_substi$Length


boxplot(positive_substi[,4:7], notch=T, outline=F, main="Positive")
boxplot(negative_substi[,4:7], notch=T, outline=F,main="Negative")
boxplot(downstream_substi[,4:7], notch=T, outline=F,main="Downstream")
boxplot(upstream_substi[,4:7], notch=T, outline=F, main="Upstream")


