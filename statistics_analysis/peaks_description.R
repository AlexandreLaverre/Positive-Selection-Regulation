library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(ROCR)
library(qvalue)
library(data.table)

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 

TF_col <- c("red", "#4DAF4A", "magenta", "black", "orange")
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"

species <- c("dog", "human", "macaca", "mouse", "rat", "spretus", "caroli", "chicken", "cat","rabbit")#cattle
TFs <- c("CTCF", "HNF6", "CEBPA", "FOXA1", "HNF4A")
divergence <- c(0.005, 6.4, 3, 2.5, 12, 2.5, 4.5, 33, 4.6, 18) #, 0.5)
names(divergence) <- species

peaks <- read.table(paste0(path, "peaks_calling/peaks_numbers.txt"),h=T, fill=T)
peaks$species <- factor(peaks$species, levels = species)
peaks$TF <- factor(peaks$TF, levels = TFs)

# Proportion signif peaks
peaks$tested_peaks <- NA
peaks$signif_low <- NA
peaks$signif_high <- NA
peaks$signif_FDR <- NA
peaks$signif_qval <- NA
CEX=1.3

#species= c("human")
for (sp in species){
  #pdf(paste0(path, "/figures/", sp, "_deltaSVM_summary.pdf"), width=8)
  for (TF in TFs){
    file = paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt")
    
    if (file.exists(file)){
      test <- read.table(file,h=T)     
      nb_signif_high <- nrow(test[which(test$pval.high < 0.05),])
      nb_signif_low <- nrow(test[which(test$pval.high > 0.95),])
      peaks[which(peaks$species==sp & peaks$TF == TF),]$tested_peaks = nrow(test)
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_high = nb_signif_high
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_low = nb_signif_low
      
      test$FDR <- p.adjust(test$pval.high, method="fdr")
      nb_signif_FDR <- nrow(test[which(test$FDR <= 0.1),])
      test$qval <- qvalue(test$pval.high)$qvalues
      nb_signif_qval <- nrow(test[which(test$FDR < 0.1),])#nrow(test[which(test$qval < 0.05 | test$qval > 0.95),])
      
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_FDR=nb_signif_FDR
      peaks[which(peaks$species==sp & peaks$TF == TF),]$signif_qval=nb_signif_qval
      
      if (sp %in% c("human", "mouse")){
        par(mfrow=c(2,2))
        hist(test$deltaSVM, breaks = 100, main=TF, xlab="deltaSVM", xlim=c(-30,30),
             cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[4])
        #hist(test$med.deltaSVM.simul, breaks = 200, main="", xlab="simul deltaSVM (median)", xlim=c(-10,10),
        #     cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[5])
        hist(test$FDR, breaks = 60, main="", xlab="qvalues", 
             cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[3])
        
        barplot(c(nb_signif_high*100/nrow(test),nb_signif_FDR*100/nrow(test), nb_signif_qval*100/nrow(test)), col=c(pal[1],pal[2], pal[3]),
                cex.lab=CEX, cex.main=CEX, cex.axis=CEX, cex.names = CEX, las=3,
                ylab="% of positive peaks", main="", names=c("pval\n0.05", "FDR\n0.1", "qval\n0.05"))
        
        ##Model performance 
        model = paste0(path, "positive_selection/", sp, "/new_run/", TF, ".cvpred.txt")
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
  }
  #dev.off()
}

peaks$signif_peaks <- peaks$signif_low + peaks$signif_high

Jialin_peaks <- read.table("/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/human/deltaSVM/CEBPA/hsap_CEBPA_deltaSVM_highertailTest.txt", h=F)    
colnames(Jialin_peaks) <- c("ID", "delta", "sub", "pval")
CEBPA_tested_peaks = nrow(Jialin_peaks)
CEBPA_signif_peaks = nrow(Jialin_peaks[which(Jialin_peaks$pval < 0.01),])
Jialin_peaks <- read.table("/Users/alaverre/Documents/Detecting_positive_selection/Tools/JialinTool/data/human/deltaSVM/HNF4A/hsap_HNF4A_deltaSVM_highertailTest.txt", h=F)    
colnames(Jialin_peaks) <- c("ID", "delta", "sub", "pval")
HNF4A_tested_peaks = nrow(Jialin_peaks)
HNF4A_signif_peaks = nrow(Jialin_peaks[which(Jialin_peaks$pval < 0.01),])


CEBPA_tested_peaks_ME = peaks[which(peaks$species == "human" & peaks$TF == "CEBPA"),]$tested_peaks
CEBPA_signif_peaks_ME = peaks[which(peaks$species == "human" & peaks$TF == "CEBPA"),]$signif_FDR
HNF4A_tested_peaks_ME = peaks[which(peaks$species == "human" & peaks$TF == "HNF4A"),]$tested_peaks
HNF4A_signif_peaks_ME = peaks[which(peaks$species == "human" & peaks$TF == "HNF4A"),]$signif_FDR

tested <- c(CEBPA_tested_peaks, HNF4A_tested_peaks, CEBPA_tested_peaks_ME, HNF4A_tested_peaks_ME)
signif <- c(CEBPA_signif_peaks, HNF4A_signif_peaks, CEBPA_signif_peaks_ME, HNF4A_signif_peaks_ME)
num_samples <- length(tested)

# Create a matrix with the values for each sample
values <- rbind(signif, tested)

# Define colors for positive and negative bars
colors <- c("green", "red")

# Create the barplot
barplot(values, beside = F, col = colors, ylim = range(signif, tested),
        names.arg = 1:num_samples, xlab = "Samples", ylab = "Count", main = "Barplot")

a <- barplot(peaks[which(peaks$species %in% c("human", "mouse")),]$signif_FDR, col=c(rep("firebrick",5), rep("lightblue", 5)),
        ylab="Number of positive peaks")

TF_names <- peaks[which(peaks$species %in% c("human", "mouse")),]$TF
text(x = a[,1], y = par("usr")[3] - 10, labels = TF_names, srt = 45, adj = 1, xpd = TRUE, cex=0.8)

legend("topleft", legend = c("Human", "Mouse"), fill = c("firebrick", "lightblue"), bty="n")

################################################################################
pdf(paste0(path,"figures/ChIP-seq_peaks_global_description.pdf"),height=6, width=7 )

# Describing samples
ggplot(peaks, aes(species, TF, fill=total_peaks)) + 
  geom_tile(colour = "black") + 
  #scale_fill_distiller(palette = "RdPu", direction=1) +
  #scale_fill_viridis(name="Nb peaks", alpha=0.9, discrete = F) +
  scale_fill_gradient(low="papayawhip", high="red", name = "# TFBS") +
  coord_equal(expand = 0) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        legend.text = element_text(size=15), legend.title = element_text(size=19),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=19), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Transcription factors", x = "Species")

# Total peaks
par(mfrow=(c(1,2)))
par(mai=c(1,0.8,0.8,0.3))
boxplot(peaks$total_peaks~peaks$species, ylab="Nb detected peaks", xlab="", las=3)
boxplot(peaks$total_peaks~peaks$TF, ylab="Nb detected peaks", xlab="", las=3)

# Proportion tested peaks
par(mfrow=(c(1,2)))
par(mai=c(1,0.8,0.4,0))
tested <- boxplot(peaks$tested_peaks/peaks$total~peaks$species, ylab="Nb tested peaks / total", xlab="", las=2)
boxplot(peaks$tested_peaks/peaks$total~peaks$TF, ylab="Nb tested peaks / total", xlab="", las=2)

col_sp = pal[1:length(divergence)] 
names(col_sp) = species
plot((100*tested$stats[3,])~divergence, col=pal[1:length(divergence)], pch=16, mgp=c(3,1,0), las=1,
     ylab="% tested peaks", xlab="Divergence time to closest species (MY)", cex.lab=1.2, cex.axis=1.2)

cor = cor.test(tested$stats[3,], divergence)
text(20, 0.2, labels=paste0("R2=", round(cor$estimate, 2)), cex=0.8)
plot.new()
legend("left", species, col=pal[1:length(divergence)], pch=16, bty='n', ncol=2, inset=c(-0.15, 0), xpd=T, cex=1.2)

par(mfrow=(c(2,2)))
par(mai=c(0.9,0.8,0.1,0.1))
boxplot(peaks$signif_peaks*100/peaks$tested_peaks~peaks$species, ylab="Proportion peaks pval<0.05", xlab="", las=2)
boxplot(peaks$signif_peaks*100/peaks$tested_peaks~peaks$TF, ylab="Proportion peaks pval<0.05", xlab="", las=2)

boxplot(peaks$signif_FDR*100/peaks$tested_peaks~peaks$species, ylab="% positive TFBS", xlab="",
        las=2, col=col_sp, medlwd=0.7, cex.axis=1.4, cex.lab=1.4)
boxplot(peaks$signif_FDR*100/peaks$tested_peaks~peaks$TF, ylab="% positive TFBS", xlab="",
        las=2, col=TF_col, medlwd=0.7, cex.axis=1.4, cex.lab=1.4,)


boxplot(peaks$signif_high*100/peaks$signif_peaks~peaks$species, ylab="Proportion positive deltaSVM", xlab="", las=2)
boxplot(peaks$signif_high*100/peaks$signif_peaks~peaks$TF, ylab="Proportion positive deltaSVM", xlab="", las=2)


#species <- c("dog", "macaca", "mouse", "cat", "human", "rat", "rabbit", "chicken")
peaks$species <- factor(peaks$species, levels = species)

tested <- boxplot((peaks$signif_qval*100)/peaks$tested_peaks~peaks$species, ylab="% positive peaks", xlab="", 
                  las=2, cex.lab=1.2, cex.axis=1.2, col=col_sp[species])

boxplot((peaks$signif_qval*100)/peaks$tested_peaks~peaks$TF, ylab="% positive peaks", xlab="", las=2, cex.lab=1.2, cex.axis=1.2)

plot(tested$stats[3,]*100~divergence, col=pal[1:length(divergence)], pch=16, mgp=c(3,1,0),
     ylab="% positive peaks", xlab="Divergence time (MY)", las=1)

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
