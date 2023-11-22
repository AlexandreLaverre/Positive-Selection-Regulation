
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
TF = "CEBPA"
species <- c("Mmusculus", "Mspretus", "Mcaroli")

################################################################################
# Datas
delta <- list()
nb_signif <- c()
for (sp in species){
  delta_SVM <- read.table(paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt"), header=T, row.names = 1)
  
  conserved <- read.table(paste0(path, "peaks_overlap/peaks/", sp, "_conserved_peaks.txt"))$V1
  gain <- read.table(paste0(path, "peaks_overlap/peaks/", sp, "_gain_peaks.txt"))$V1
  loss <- read.table(paste0(path, "peaks_overlap/peaks/", sp, "_loss_peaks.txt"))$V1
  
  conserved <- conserved[conserved %in% row.names(delta_SVM)]
  gain <- gain[gain %in% row.names(delta_SVM)]
  loss <- loss[loss %in% row.names(delta_SVM)]
  
  delta_SVM_gain <- delta_SVM[gain,]
  delta_SVM_conserved <- delta_SVM[conserved,]
  
  nb_gain_signif = nrow(delta_SVM_gain[which(delta_SVM_gain$pval.high<0.05),])
  nb_conserv_signif = nrow(delta_SVM_conserved[which(delta_SVM_conserved$pval.high<0.05),])
  
  delta[[paste(sp, "gain")]] <- delta_SVM_gain$deltaSVM
  delta[[paste(sp, "conserved")]] <- delta_SVM_conserved$deltaSVM
  delta[[paste(sp, "NA")]] <- NA
  
  nb_signif <- c(nb_signif, nb_gain_signif, nb_conserv_signif, NA)
}

################################################################################
# Plots
type = c("Gain", "Conserved")
col = c("lawngreen", "skyblue3", NA)

pdf(file=paste0(path, "figures/deltaSVM_mice_gain_vs_conserved.pdf"),height=6, width=8)
par(mfrow=c(1,2), mgp=c(2.3,1,0))
par(mai=c(0.8,0.75,0.8,0.3), xpd=T)
a <- boxplot(delta[-length(delta)], ylab="deltaSVM", names=F, col=col, outline=F, notch=T, las=1, cex.lab=1.2)
mtext(species, side=1, at=c(1.5, 4.5, 7.5), cex=1.2, line=1)

barplot(nb_signif, ylab="Nb positively selected", names="", col=col, cex.lab=1.2)
legend("topright", legend=type, fill=col, bty="n", cex=1.1, inset=c(-0.1,0), x.intersp = 0.3)
mtext(species, side=1, at=c(1.3, 4.9, 8.5), cex=1.2, line=1)

dev.off()
################################################################################
