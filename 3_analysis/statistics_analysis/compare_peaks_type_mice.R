path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
TFs = c("CEBPA", "HNF4A", "FOXA1")
species <- c("Mmusculus", "Mspretus") #, "Mcaroli")
col = c("forestgreen", "orange", "deepskyblue4")

treshold = 0.2
minMut = 4

pdf(file=paste0(path, "../../results/figures/mice_gain_loss_conserved_FDR_", treshold, "_", minMut, "mut.pdf"),height=6, width=8)
par(mfrow=c(2,2), mgp=c(2.3,1,0))
par(mfrow=c(1,2))
for (TF in TFs){
  message(TF)
  par(mfrow=c(2,2))
  ################################################################################
  # Datas
  delta <- list()
  nb_signif <- c()
  nb_MLL <- c()
  nb_MLL_absolute <- c()
  
  for (sp in species){
    common_name=ifelse(sp=="Mmusculus", "mouse", ifelse(sp=="Mspretus", "spretus", "caroli"))
    sample=ifelse(sp=="Mmusculus", "Wilson", "Stefflova")
    
    delta_SVM <- read.table(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/", sample, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt"), header=T, row.names = 1)
    delta_SVM$pval.two.tailed <- apply(delta_SVM, 1, function(x) 2*min(as.numeric(x["pval.high"]), 1-as.numeric(x["pval.high"])))
    delta_SVM$pval.high <- delta_SVM$pval.two.tailed
    delta_SVM <- delta_SVM[which(delta_SVM$NbSub >= minMut),]
    delta_SVM$FDR <- p.adjust(delta_SVM$pval.high, method="fdr")
    
    delta_loss <- read.table(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/loss/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt"), header=T, row.names = 1)
    delta_loss$pval.two.tailed <- apply(delta_loss, 1, function(x) 2*min(as.numeric(x["pval.high"]), 1-as.numeric(x["pval.high"])))
    delta_loss$pval.high <- delta_loss$pval.two.tailed
    delta_loss <- delta_loss[which(delta_loss$NbSub >= minMut),]
    delta_loss$FDR <- p.adjust(delta_loss$pval.high, method="fdr")
    
    conserved <- read.table(paste0(path, "peaks_overlap/NarrowPeaks/mouse/", TF, "/peaks_by_class/", sp, "_conserved_peaks.txt"))$V4
    gain <- read.table(paste0(path, "peaks_overlap/NarrowPeaks/mouse/", TF, "/peaks_by_class/", sp, "_gain_peaks.txt"))$V4
    loss <- read.table(paste0(path, "peaks_overlap/NarrowPeaks/mouse/", TF, "/peaks_by_class/", sp, "_loss_peaks.txt"))$V4
    
    conserved <- conserved[conserved %in% row.names(delta_SVM)]
    gain <- gain[gain %in% row.names(delta_SVM)]
    loss <- loss[loss %in% row.names(delta_loss)]
    
    delta_SVM_gain <- delta_SVM[gain,]
    delta_SVM_conserved <- delta_SVM[conserved,]
    delta_SVM_loss <- delta_loss[loss,]
    
    nb_gain_signif = nrow(delta_SVM_gain[which(delta_SVM_gain$FDR<treshold),])/nrow(delta_SVM_gain)
    nb_conserv_signif = nrow(delta_SVM_conserved[which(delta_SVM_conserved$FDR<treshold),])/nrow(delta_SVM_conserved)
    nb_loss_signif = nrow(delta_SVM_loss[which(delta_SVM_loss$FDR<treshold),])/nrow(delta_SVM_loss)
    
    delta[[paste(sp, "gain")]] <- delta_SVM_gain$deltaSVM
    delta[[paste(sp, "conserved")]] <- delta_SVM_conserved$deltaSVM
    delta[[paste(sp, "loss")]] <- delta_SVM_loss$deltaSVM
    nb_signif <- c(nb_signif, nb_gain_signif, nb_conserv_signif, nb_loss_signif, NA)
    
    # MLL ranked
    MLL <- read.csv(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/", sample, "/", TF, "/Tests/MLE_summary_exact.csv"), header=T, row.names = 1)
    MLL_loss <- read.csv(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/loss/", TF, "/Tests/MLE_summary_exact_ranked.csv"), header=T, row.names = 1)
    MLL <- MLL[which(MLL$Nmut >= minMut),]
    MLL_loss <- MLL_loss[which(MLL_loss$Nmut >= minMut),]
    
    MLL$FDR <- p.adjust(MLL$p_value_purif_pos, method="fdr")
    MLL_loss$FDR <- p.adjust(MLL_loss$p_value_purif_pos, method="fdr")
    
    MLL$Conclusion <- ifelse(MLL$FDR > treshold , "Neutral model", MLL$Conclusion)
    MLL_loss$Conclusion <- ifelse(MLL_loss$FDR > treshold , "Neutral model", MLL_loss$Conclusion)
    
    conserved <- conserved[conserved %in% row.names(MLL)]
    gain <- gain[gain %in% row.names(MLL)]
    loss <- loss[loss %in% row.names(MLL_loss)]
    
    MLL_gain <- MLL[gain,]
    MLL_conserved <- MLL[conserved,]
    MLL_loss <- MLL_loss[loss,]
    
    nb_gain_signif = nrow(MLL_gain[which(grepl("Directional", MLL_gain$Conclusion)),])/nrow(MLL_gain)
    nb_conserv_signif = nrow(MLL_conserved[which(grepl("Directional", MLL_conserved$Conclusion)),])/nrow(MLL_conserved)
    nb_loss_signif = nrow(MLL_loss[which(grepl("Directional", MLL_loss$Conclusion)),])/nrow(MLL_loss)
    nb_MLL <- c(nb_MLL, nb_gain_signif, nb_conserv_signif, nb_loss_signif, NA)
    
    # MLL absolute
    MLL <- read.csv(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/", sample, "/", TF, "/Tests/MLE_summary_exact_absolute.csv"), header=T, row.names = 1)
    MLL_loss <- read.csv(paste0(path, "positive_selection/NarrowPeaks/", common_name, "/loss/", TF, "/Tests/MLE_summary_exact_absolute.csv"), header=T, row.names = 1)
    
    MLL <- MLL[which(MLL$Nmut >= minMut),]
    MLL_loss <- MLL_loss[which(MLL_loss$Nmut >= minMut),]
    
    MLL$FDR <- p.adjust(MLL$p_value_purif_pos, method="fdr")
    MLL_loss$FDR <- p.adjust(MLL_loss$p_value_purif_pos, method="fdr")
    
    MLL$Conclusion <- ifelse(MLL$FDR > treshold , "Neutral model", MLL$Conclusion)
    MLL_loss$Conclusion <- ifelse(MLL_loss$FDR > treshold , "Neutral model", MLL_loss$Conclusion)
    
    conserved <- conserved[conserved %in% row.names(MLL)]
    gain <- gain[gain %in% row.names(MLL)]
    loss <- loss[loss %in% row.names(MLL_loss)]
    
    MLL_gain <- MLL[gain,]
    MLL_conserved <- MLL[conserved,]
    MLL_loss <- MLL_loss[loss,]
    
    nb_gain_signif = nrow(MLL_gain[which(grepl("Directional", MLL_gain$Conclusion)),])/nrow(MLL_gain)
    nb_conserv_signif = nrow(MLL_conserved[which(grepl("Directional", MLL_conserved$Conclusion)),])/nrow(MLL_conserved)
    nb_loss_signif = nrow(MLL_loss[which(grepl("Directional", MLL_loss$Conclusion)),])/nrow(MLL_loss)
    nb_MLL_absolute <- c(nb_MLL_absolute, nb_gain_signif, nb_conserv_signif, nb_loss_signif, NA)
    
    a <- prop.table(table(MLL_loss[which(grepl("Directional", MLL_loss$Conclusion)),]$Conclusion))
    message(paste(signif(a[1]*100, digits=3), "% of loss peaks are Directional (-)"))
    
    MLL_final <- rbind(MLL, MLL_loss)
    MLL_final$status <- "null"
    MLL_final[gain,]$status <- "Gain"
    MLL_final[conserved,]$status <- "Conserved"
    MLL_final[loss,]$status <- "Loss"
    
    MLL_final$Conclusion <- factor(MLL_final$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))
    Conclu.status <- prop.table(table(MLL_final$Conclusion, MLL_final$status), margin = 2)
    barplot(Conclu.status[1:3,1:3]*100, col=col, las=1, ylab="Proportion MLL (%)", xlab="", main=paste(sp, TF), cex.lab=1, cex.names = 1.2)
    
    if (sp == "Mspretus"){legend("topleft", legend=c("Directional (+)", "Directional (-)", "Stabilising"), fill=col, bty="n")}
    
    # Permut
    delta_SVM_final <- rbind(delta_SVM, delta_loss)
    delta_SVM_final$status <- "null"
    delta_SVM_final[gain,]$status <- "Gain"
    delta_SVM_final[conserved,]$status <- "Conserved"
    delta_SVM_final[loss,]$status <- "Loss"
    delta_SVM_final$signif <- as.factor(ifelse(delta_SVM_final$pval.high > 0.01, "Neutral", ifelse(delta_SVM_final$deltaSVM >0, "Directional (+)", "Directional (-)")))
    delta_SVM_final$signif <- factor(delta_SVM_final$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral", "Disruptive"))
    
    Conclu.status <- prop.table(table(delta_SVM_final$signif, delta_SVM_final$status), margin = 2)
    barplot(Conclu.status[1:3,1:3]*100, col=col, las=1, ylab="Proportion Permut (%)", xlab="", main=paste(sp, TF), cex.lab=1.2, cex.names = 1.2)
    
  }
  
  ################################################################################
  # Plots
  #type = c("Gain", "Conserved", "Loss")
  #col = c("forestgreen", "orange", "deepskyblue4")
  par(mfrow=c(1,1))
  #par(mai=c(0.8,0.75,0.8,0.3), xpd=F)
  a <- boxplot(delta, ylab="Phenotypic change", names=F, col=col, outline=F, notch=T, las=1, cex.lab=1.2, main=TF)
  legend("bottomleft", legend=type, fill=col, bty="n", cex=1.1)
  mtext(species, side=1, at=c(2, 5, 8), cex=1, line=1)
  
  #barplot(nb_signif, ylab="Prop. Permut", names="", col=c(col,"white"), cex.lab=1.2)
  #legend("bottomleft", legend=type, fill=col, bty="n", cex=1.1, inset=c(-0.1,0), x.intersp = 0.3)
  #mtext(species, side=1, at=c(2, 7, 11.5), cex=1, line=1)
  
  #barplot(nb_MLL, ylab="Prop. MLL ranked", names="", col=c(col,"white"), cex.lab=1.2, ylim=c(0, 0.04))
  #mtext(species, side=1, at=c(2, 7, 11.5), cex=1, line=1)
  
  #barplot(nb_MLL_absolute, ylab="Prop. MLL absolute", names="", col=c(col,"white"), cex.lab=1.2, ylim=c(0, 0.04))
  #mtext(species, side=1, at=c(2, 7, 11.5), cex=1, line=1)

  ################################################################################
}

dev.off()
