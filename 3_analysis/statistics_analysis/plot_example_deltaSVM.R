path="/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
sp="human"
samp="Wilson"
TF="HNF4A"

# When simul
simul=TRUE
bin = "exact_ranked_50"
data = "beta_0.1Null_25Stab_25Pos"
sel = "positive"

ref=ifelse(simul, "focal", "ancestral")
obsdelta=ifelse(simul, paste0("simul_", data, "_", sel), "ancestral_to")
method=paste0("simulated_", data, "_", sel, "_", bin, "bins_threshold_0.01")

# MLE
MLE.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/MLE_summary_", method, ".csv2")
MLE <- read.csv(MLE.file, h=T, row.names = 1)
rownames(MLE) <- gsub("_.*:\\d+:\\d+:", "_", rownames(MLE))

# Observed deltas
delta.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/deltas/", obsdelta, "_observed_deltaSVM.txt")
obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col) 
row.names(deltas) <- deltas$seq_name
rownames(deltas) <- gsub("_.*:\\d+:\\d+:", "_", rownames(deltas))
deltas <- deltas[row.names(MLE),]

# All deltas 
Alldelta.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/deltas/", ref, "_all_possible_deltaSVM.txt")
Alldeltas <- read.table(Alldelta.file, h=T, sep="\t", quote="", fill=T, row.names = 1)
rownames(Alldeltas) <- gsub("_.*:\\d+:\\d+:", "_", rownames(Alldeltas))

# Distrib pval
par(mfrow=c(2,2))
hist(MLE$p_value_null_purif, breaks=100, main="Null vs Stabilising", xlab="p-values")
hist(MLE$p_value_null_pos, breaks=100, main="Null vs Directional", xlab="p-values")
hist(MLE$p_value_purif_pos, breaks=100, main="Stabilising vs Directional", xlab="p-values")


# FDR
MLE$FDR_null_pos <- p.adjust(MLE$p_value_null_pos, method="fdr")
MLE$FDR_null_purif <- p.adjust(MLE$p_value_null_purif, method="fdr")
MLE$FDR_purif_pos <- p.adjust(MLE$p_value_purif_pos, method="fdr")

#IDs = row.names(stab[which(stab$Nmut==20),])[1:4]
#special_case = rownames(special_case[which(special_case$species=="human" & special_case$TF == "CEBPA"),])
#IDs = gsub("human_CEBPA.", "", special_case)
IDs = c("chr1:32760223:32760419_Interval_1065", "chr15:70078718:70079076_Interval_19597",
        "chr1:204260819:204261152_Interval_4824", "chr1:62869418:62869692_Interval_2060")
#IDs = diff[which(diff$species=="human" & diff$TF=="CEBPA"),]$ID[1:24]
par(xpd = F)
par(mfrow=c(2,2))
for (ID in IDs){
  #ID="chr2:31837530:31837664"
  FDR=MLE[ID,]$FDR_null_pos
  pval=MLE[ID,]$p_value_null_pos
  hist(as.numeric(Alldeltas[ID,]), breaks=200,  
       main=paste("FDR=",signif(FDR, 2), "pval=", signif(pval, 2)), xlab="deltaSVM")
  abline(v=as.numeric(deltas[ID,5:100]), col="red")
}



##### Plot Beta distribution 
alpha=5
beta=5


### 
p = seq(0, 1, length=1000)
#create plot of Beta distribution with shape parameters 2 and 10
plot(p, dbeta(p, alpha, beta), type='l', ylab="Relative Fixation Probability", xlab=expression("scaled "*Delta*"SVM"), 
     las=1, col="forestgreen", lwd=3, cex.axis=1, cex.lab=1.1, main="Beta distribution", axes=T)
legend("topleft", legend=expression(alpha*"="*beta*"=5"), bty="n", cex=1.3)

lines(p, dbeta(p, 40, 2), col="navy", lwd=3)
lines(p, rep(5, length(p)), col="firebrick", lwd=3)
axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(-5, -2.5, 0, 2.5, 5))
axis(2, at=c(0, 2.5, 5, 7.5, 10), labels=c(0, 0.5, 1, 1.5, 2), las=2)
legend("topleft", legend=c("Beta=2", "Beta=1"), lty=c(1,1,1), lwd=3,
       col=c("navy", "forestgreen"), cex=1.2)



