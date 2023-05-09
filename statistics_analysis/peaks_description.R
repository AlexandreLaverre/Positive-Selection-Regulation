library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(ROCR)
library(qvalue)
library(data.table)

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"

species <- c("dog", "human", "macaca", "mouse", "rat", "chicken", "cat","rabbit", "cattle")
TFs <- c("CEBPA", "FOXA1", "HNF4A", "HNF6")
divergence <- c(0.005, 6.4, 3, 2.5, 12, 33, 4.6, 18, 0.5)
names(divergence) <- species

peaks <- read.table(paste0(path, "peaks_calling/peaks_numbers.txt"),h=T, fill=T)
peaks$species <- factor(peaks$species, levels = species)
peaks$TF <- factor(peaks$TF, levels = TFs)

# Proportion signif peaks
peaks$tested_peaks <- NA
peaks$signif_low <- NA
peaks$signif_high <- NA
peaks$signif_FDR <- NA
CEX=1.3

species= c("human")
for (sp in species){
  pdf(paste0(path, "/figures/", sp, "_deltaSVM_summary.pdf"), width=8)
  for (TF in TFs){
    file = paste0(path, "positive_selection/", sp, "/new_run/", TF, "_PosSelTest_deltaSVM_10000permutations.txt")
    model = paste0(path, "positive_selection/", sp, "/new_run/", TF, ".cvpred.txt")
    if (file.exists(file)){
      test <- read.table(file,h=T)     
      nb_signif_high <- nrow(test[which(test$pval.high < 0.05),])
      nb_signif_low <- nrow(test[which(test$pval.high > 0.95),])
      peaks[which(peaks$species==sp & peaks$TF == TF),]$tested_peaks = nrow(test)
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_high = nb_signif_high
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_low = nb_signif_low
      
      test$FDR <- p.adjust(test$pval.high, method="fdr")
      nb_signif_FDR <- nrow(test[which(test$FDR < 0.1),])
      test$qval <- qvalue(test$pval.high)$qvalues
      nb_signif_qval <- nrow(test[which(test$qval < 0.05),])
      
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_FDR=nb_signif_FDR
      
      par(mfrow=c(2,2))
      hist(test$deltaSVM, breaks = 100, main=TF, xlab="deltaSVM", xlim=c(-30,30),
           cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[4])
      hist(test$med.deltaSVM.simul, breaks = 200, main="", xlab="simul deltaSVM (median)", xlim=c(-10,10),
           cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[5])
      hist(test$qval, breaks = 60, main="", xlab="qvalues", 
           cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[3])
      
      barplot(c(nb_signif_high*100/nrow(test),nb_signif_FDR*100/nrow(test), nb_signif_qval*100/nrow(test)), col=c(pal[1],pal[2], pal[3]),
              cex.lab=CEX, cex.main=CEX, cex.axis=CEX, cex.names = CEX, las=3,
              ylab="% of positive peaks", main="", names=c("pval\n0.05", "FDR\n0.1", "qval\n0.05"))
      
      ##Model performance 
      cv <- fread(model)
      colnames(cv)<-c("position","prediction","real_state","cv_number")
      
      pred <- prediction(cv$prediction, cv$real_state) 
      perf <- performance(pred, "tpr", "fpr" )
      
      par(mfrow=c(1,1))
      plot(perf,lwd = 3,cex=1.4, main=paste(sp, TF, "prediction"))
      
      auc_result <-performance(pred, measure = "auc")
      AUC <- signif(unlist(slot(auc_result, "y.values")), 3)
      legend(0.3,0.6, paste0("(AUC = ",  AUC, ")"), bty="n",lwd = 3,cex=1.4) 

    }
  }
  dev.off()
}

peaks$signif_peaks <- peaks$signif_low + peaks$signif_high

################################################################################
pdf(paste0(path,"figures/ChIP-seq_peaks_global_description.pdf"),height=6, width=7 )

# Describing samples
ggplot(peaks, aes(species, TF, fill=total_peaks)) + 
  geom_tile(colour = "black") + 
  #scale_fill_distiller(palette = "RdPu", direction=1) +
  #scale_fill_viridis(name="Nb peaks", alpha=0.9, discrete = F) +
  scale_fill_gradient(low="papayawhip", high="red", name = "Nb peaks") +
  coord_equal(expand = 0) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=12), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Transcription factors", x = "Species")

# Total peaks
par(mfrow=(c(1,2)))
par(mai=c(1,0.8,0.8,0.3))
boxplot(peaks$total_peaks~peaks$species, ylab="Nb detected peaks", xlab="", las=3)
boxplot(peaks$total_peaks~peaks$TF, ylab="Nb detected peaks", xlab="", las=3)

# Proportion tested peaks
par(mfrow=(c(2,2)))
par(mai=c(0.6,0.8,0.4,0.3))
tested <- boxplot(peaks$tested_peaks/peaks$total~peaks$species, ylab="Nb tested peaks / total", xlab="", las=2)
boxplot(peaks$tested_peaks/peaks$total~peaks$TF, ylab="Nb tested peaks / total", xlab="", las=2)
plot(tested$stats[3,]~divergence, col=pal[1:length(divergence)], pch=16, mgp=c(2,1,0),
     ylab="Nb tested peaks/total (median)", xlab="Divergence time (MY)")

cor = cor.test(tested$stats[3,], divergence)
text(25, 0.2, labels=paste0("R2=", round(cor$estimate, 2)), cex=0.8)
plot.new()
legend("left", species, col=pal[1:length(divergence)], pch=16, bty='n', ncol=2, inset=c(-0.35, 0), xpd=T)

par(mfrow=(c(1,2)))
par(mai=c(1,0.8,0.8,0.3))
boxplot(peaks$signif_peaks/peaks$tested_peaks~peaks$species, ylab="Proportion peaks pval<0.05", xlab="", las=2)
boxplot(peaks$signif_peaks/peaks$tested_peaks~peaks$TF, ylab="Proportion peaks pval<0.05", xlab="", las=2)

boxplot(peaks$signif_FDR/peaks$tested_peaks~peaks$species, ylab="Proportion peaks FDR<0.1", xlab="", las=2)
boxplot(peaks$signif_FDR/peaks$tested_peaks~peaks$TF, ylab="Proportion peaks FDR<0.1", xlab="", las=2)

boxplot(peaks$signif_high/peaks$signif_peaks~peaks$species, ylab="Proportion positive deltaSVM", xlab="", las=2)
boxplot(peaks$signif_high/peaks$signif_peaks~peaks$TF, ylab="Proportion positive deltaSVM", xlab="", las=2)

dev.off()

################################################################################
species=c("mouse", "human")
TFs=c("CEBPA", "HNF4A", "HNF6", "FOXA1", "CTCF")

for (sp in species){
  pdf(paste0(path, "/figures/", sp, "_deltaSVM_phyloP.pdf"), width=8)
  suffix = ifelse(sp == "human", "_17way.txt", "_60way.glire.txt")
  for (TF in TFs){
    
    PosSel <- read.table(paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt"),h=T)
    phyloP <- read.table(paste0(path, "phyloP/", sp, "/", TF, suffix), h=T)
    rownames(PosSel) <- PosSel$ID
    rownames(phyloP) <- phyloP$ID
    phyloP <- phyloP[rownames(PosSel),]
    
    phyloP$class <- cut(phyloP$Score, quantile(phyloP$Score, probs = seq(0, 1, 0.2), na.rm=T), include.lowest=T, 
                        labels=c("20%", "40%", "60%", "80%", "100%"))
    
    par(mfrow=c(2,2))
    #deltaSVM
    plot(PosSel$deltaSVM~phyloP$Score, cex=0.1, main=paste(sp, TF), xlab="phyloP score", ylab="deltaSVM")
    cor = cor.test(phyloP$Score, PosSel$deltaSVM)
    mtext(paste0("R2=", round(cor$estimate, 2)), cex=0.8)
    
    boxplot(PosSel$deltaSVM~phyloP$class, notch=T, xlab="phyloP score quantile", ylab="deltaSVM", las=1, outline=F)
    
    #deltaSVM pval
    plot(PosSel$pval.high~phyloP$Score, cex=0.1, main="", xlab="phyloP score", ylab="deltaSVM pvalue")
    cor = cor.test(phyloP$Score, PosSel$pval.high)
    
    mtext(paste0("R2=", round(cor$estimate, 2)), cex=0.8)
    boxplot(PosSel$pval.high~phyloP$class, notch=T, xlab="phyloP score quantile", ylab="deltaSVM pvalue", las=1)
  }
  dev.off()
}

################################################################################
