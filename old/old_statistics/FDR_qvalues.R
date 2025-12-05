sp="caroli"
TF="CEBPA"
data = all_MLE[which(all_MLE$species==sp & all_MLE$TF==TF),]

pdf(paste0(pathFigures, sp, "_", TF, "_pvalues_analyses.pdf"))
    
# Distribution pvalues vs qvalues
par(mfrow=c(3,2))
hist(data$pval.high, breaks=50, main=paste("Permut one-tailed"), xlab="pval. Permut")
qval= qvalue(data$pval.high)$qvalues
hist(qval, breaks=20, main="", xlab="qval. Permut", xlim=c(0,1))

hist(data$pval.two.tailed, breaks=50, main=paste("Permut two-tailed"), xlab="pval. Permut")
qval= qvalue(data$pval.two.tailed)$qvalues
hist(qval, breaks=50, main="", xlab="qval. Permut", xlim=c(0,1))

hist(data$p_value_null_pos, breaks=100, main=paste("Test MLL"), xlab="pval. MLL null_pos")
qval= qvalue(data$p_value_null_pos)$qvalues
hist(qval, breaks=50, main="", xlab="qval. MLL null_pos", xlim=c(0,1))

# Correlation between Tests
data$class_pval <- cut(data$p_value_null_pos, breaks=quantile(data$p_value_null_pos, probs = seq(0, 1, 0.05)), include.lowest = T)
data$class_pval_two <- cut(data$pval.two.tailed, breaks=quantile(data$pval.two.tailed, probs = seq(0, 1, 0.05)), include.lowest = T)
data$class_pval_high <- cut(data$pval.high, breaks=quantile(data$pval.high, probs = seq(0, 1, 0.05)), include.lowest = T)

par(mfrow=c(2,2))
boxplot(data$p_value_null_pos~data$class_pval_high, outline=F, notch=T, xlab="Quantile of pval one-tailed", ylab="MLL pval Null vs Dir")
legend("topleft", legend=paste("R=",round(cor.test(data$p_value_null_pos, data$pval.high)$estimate, digits=3)), bty="n")
boxplot(data$p_value_null_pos~data$class_pval_two, outline=F, notch=T, xlab="Quantile of pval two-tailed", ylab="MLL pval Null vs Dir")
legend("topleft", legend=paste("R=",round(cor.test(data$p_value_null_pos, data$pval.two.tailed)$estimate, digits=3)), bty="n")

boxplot(data$pval.high~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="pval one-tailed")
legend("topleft", legend=paste("R=",round(cor.test(data$pval.high, data$p_value_null_pos)$estimate, digits=3)), bty="n")
boxplot(data$pval.two.tailed~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="pval two-tailed")
legend("topleft", legend=paste("R=",round(cor.test(data$pval.two.tailed, data$p_value_null_pos)$estimate, digits=3)), bty="n")

# Correlations with sequences features
par(mfrow=c(2,2))
boxplot(data$Nmut~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="Nmut", main="MLL pvalues")
legend("topright", legend=paste("R=",round(cor.test(data$Nmut, data$p_value_null_pos)$estimate, digits=3)), bty="n")
boxplot(data$deltaSVM~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="deltaSVM")
legend("topright", legend=paste("R=",round(cor.test(data$deltaSVM, data$p_value_null_pos)$estimate, digits=3)), bty="n")
boxplot(data$SVM~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="SVM")
legend("topright", legend=paste("R=",round(cor.test(data$SVM, data$p_value_null_pos)$estimate, digits=3)), bty="n")
boxplot(data$AlphaPos~data$class_pval, outline=F, notch=T, xlab="Quantile of MLL pval", ylab="AlphaPos")
legend("topright", legend=paste("R=",round(cor.test(data$AlphaPos, data$p_value_null_pos)$estimate, digits=3)), bty="n")

par(mfrow=c(2,2))
boxplot(data$Nmut~data$class_pval_two, outline=F, notch=T, xlab="Quantile of pval two-tailed", ylab="Nmut", main="Permut 2-tailed pvalues")
legend("topright", legend=paste("R=",round(cor.test(data$Nmut, data$pval.two.tailed)$estimate, digits=3)), bty="n")
boxplot(data$deltaSVM~data$class_pval_two, outline=F, notch=T, xlab="Quantile of pval two-tailed", ylab="deltaSVM")
legend("topright", legend=paste("R=",round(cor.test(data$deltaSVM, data$pval.two.tailed)$estimate, digits=3)), bty="n")
boxplot(data$SVM~data$class_pval_two, outline=F, notch=T, xlab="Quantile of pval two-tailed", ylab="SVM")
legend("topright", legend=paste("R=",round(cor.test(data$SVM, data$pval.two.tailed)$estimate, digits=3)), bty="n")
boxplot(data$AlphaPos~data$class_pval_two, outline=F, notch=T, xlab="Quantile of pval two-tailed", ylab="AlphaPos")
legend("topright", legend=paste("R=",round(cor.test(data$AlphaPos, data$pval.two.tailed)$estimate, digits=3)), bty="n")


# Proportion significant according to treshold
tresholds=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
permut <- c()
permut.qval <- c()
permut.fdr <- c()
permut.two <- c()
permut.two.qval <- c()
permut.two.fdr <- c()
LRT <- c()
LRT.qval <- c()
LRT.fdr <- c()

for (t in tresholds){
  permut.two <- c(permut.two, sum(data$pval.two.tailed < t)/nrow(data))
  permut.two.qval <- c(permut.two.qval, sum(qvalue(data$pval.two.tailed)$qvalues < t)/nrow(data))
  permut.two.fdr <- c(permut.two.fdr, sum(p.adjust(data$pval.two.tailed, method="fdr") < t)/nrow(data))
  
  permut <- c(permut, sum(data$pval.high < t)/length(data$pval.high))
  permut.qval <- c(permut.qval, sum(qvalue(data$pval.high)$qvalues < t)/nrow(data))
  permut.fdr <- c(permut.fdr, sum(p.adjust(data$pval.high, method="fdr") < t)/nrow(data))
  
  LRT <- c(LRT, sum(data$p_value_null_pos < t )/nrow(data))
  LRT.qval <- c(LRT.qval, sum(qvalue(data$p_value_null_pos)$qvalues < t)/nrow(data))
  LRT.fdr <- c(LRT.fdr, sum(p.adjust(data$p_value_null_pos, method="fdr") < t)/nrow(data))
}

par(mfrow=c(1,3))
plot(LRT, xaxt="n", type="b", ylim=c(0, max(LRT, permut, permut.two)),
     xlab="Treshold", ylab="Proportion Directional", main="P-values")
axis(1, at=1:6,labels=tresholds)
lines(permut, type="b", col="red")
lines(permut.two, type="b", col="blue")
legend("topleft", legend=c("MLL", "Permut 1-tail", "Permut 2-tail"), fill=c("black","red", "blue"), bty="n")

plot(LRT.qval, xaxt="n", type="b", ylim=c(0, max(permut.qval, LRT.qval, permut.two.qval)),
     xlab="Treshold", ylab="Proportion Directional", main="Q-values")
axis(1, at=1:6,labels=tresholds)
lines(permut.qval, type="b", col="red")
lines(permut.two.qval, type="b", col="blue")
legend("topleft", legend=c("MLL", "Permut 1-tail", "Permut 2-tail"), fill=c("black","red", "blue"), bty="n")

plot(LRT.fdr, xaxt="n", type="b", ylim=c(0, max(LRT.fdr, permut.fdr, permut.two.fdr)),
     xlab="Treshold", ylab="Proportion Directional", main="FDR")
axis(1, at=1:6,labels=tresholds)
lines(permut.fdr, type="b", col="red")
lines(permut.two.fdr, type="b", col="blue")
legend("topleft", legend=c("MLL", "Permut 1-tail", "Permut 2-tail"), fill=c("black","red", "blue"), bty="n")

dev.off()