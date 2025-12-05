library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)

path = "MLE_SNP_human/"
TFs = c("CEBPA", "FOXA1", "HNF4A", "HNF6")


for (TF in TFs){
  #pdf(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/Polymorphism_analyse_", TF, ".pdf"))
  MLE <- read.table(paste0(path, "/human_", TF, "_MLE_summarys.csv"), h=T, sep=",")
  rownames(MLE) <- MLE$ID
  
  SNP <- read.table(paste0(path, "/human_", TF, "_SNP_SelectionCoefficient.txt"), h=T, sep="\t")
  SNP <- SNP[which(SNP$NbAlt>1),]
  
  SNP$freq <- log(SNP$NbAlt/SNP$NbTot)
  SNP$classDeltaSVM <- cut(SNP$DeltaSVM, breaks=c(min(SNP$DeltaSVM), -5, -1, 0, 1, 4, max(SNP$DeltaSVM)), include.lowest = T)
  SNP$classSelCoefStab <- cut(SNP$SelCoefStab, breaks=quantile(SNP$SelCoefStab), include.lowest = T)
  SNP$classSelCoefPos <- cut(SNP$SelCoefPos, breaks=quantile(SNP$SelCoefPos), include.lowest = T)
  
  par(mfrow=c(2,2))
  names=c(paste0("Alt=Ancestral\nN=",table(SNP$Flag)[1]), paste0("Ref=Ancestral\nN=",table(SNP$Flag)[2])) 
  boxplot(SNP$freq~SNP$Flag, notch=T, ylab="Alternative Allele Frequency (log)", xlab="",
          names=names, main="Polarisation issue?")
  boxplot(SNP$DeltaSVM~SNP$Flag, notch=T, ylab="DeltaSVM", xlab="", outline=F, names=c("Alt=Ancestral","Ref=Ancestral"))
  boxplot(SNP$SelCoefStab~SNP$Flag, notch=T, ylab="Stabilising Selection Coefficient", xlab="", outline=F, names=c("Alt=Ancestral","Ref=Ancestral"))
  boxplot(SNP$SelCoefPos~SNP$Flag, notch=T, ylab="Directional Selection Coefficient", xlab="", outline=F, names=c("Alt=Ancestral","Ref=Ancestral"))
  
  SNP <- SNP[which(SNP$Flag == "ref_ancestral"),]
  
  par(mfrow=c(2,2))
  boxplot(SNP$freq~SNP$classDeltaSVM, notch=T, outline=F, xlab="deltaSVM", ylab="Alternative Allele Frequency (log)", las=2)
  boxplot(SNP$SelCoefStab~SNP$classDeltaSVM, notch=T, outline=F, xlab="deltaSVM", ylab="Stabilising Selection Coefficient", las=2)
  boxplot(SNP$SelCoefPos~SNP$classDeltaSVM, notch=T, outline=F, xlab="deltaSVM", ylab="Directional Selection Coefficient", las=2)
  
  SNP_pos <- SNP[which(SNP$ID %in% MLE[which(MLE$AlphaPos>MLE$BetaPos),"ID"]),]
  SNP_neg <- SNP[which(SNP$ID %in% MLE[which(MLE$AlphaPos<MLE$BetaPos),"ID"]),]
  
  boxplot(SNP_pos$SelCoefPos~SNP_pos$classDeltaSVM, notch=T, outline=F, xlab="deltaSVM", ylab="Directional Selection Coefficient", las=2)
  boxplot(SNP_neg$SelCoefPos~SNP_neg$classDeltaSVM, notch=T, outline=F, xlab="deltaSVM", ylab="Directional Selection Coefficient", las=2)
  
  par(mfrow=c(1,2))
  boxplot(SNP$freq~SNP$classSelCoefStab, notch=T, outline=F, ylab="Alternative Allele Frequency (log)", xlab="Stabilising Selection Coefficient")
  boxplot(SNP$freq~SNP$classSelCoefPos, notch=T, outline=F, ylab="Alternative Allele Frequency (log)", xlab="Directional Selection Coefficient")
  
  SNP$Conclusion <- MLE[SNP$ID, "Conclusion"]
  boxplot(SNP$freq~SNP$Conclusion, notch=T, outline=F)
  boxplot(SNP$SelCoefStab~SNP$Conclusion, notch=T, outline=F)
  boxplot(SNP$SelCoefPos~SNP$Conclusion, notch=T, outline=F)
  boxplot(SNP$DeltaSVM~SNP$Conclusion, notch=T, outline=F)
  
  SNP_dir_pos <- SNP_pos[which(SNP_pos$Conclusion=="Positive model"),]
  SNP_dir_neg <- SNP_neg[which(SNP_neg$Conclusion=="Positive model"),]
  SNP_neut <- SNP[which(SNP$Conclusion=="Neutral model"),]
  SNP_stab <- SNP[which(SNP$Conclusion=="Stabilizing model"),]
  subsamp <- list(SNP_dir_pos, SNP_dir_neg, SNP_neut, SNP_stab)
  
  for (samp in subsamp){
    par(mfrow=c(2,2))
    boxplot(samp$freq~samp$classDeltaSVM, notch=T, outline=F, ylab="Alternative Allele Frequency (log)", xlab="deltaSVM")
    boxplot(samp$freq~samp$classSelCoefStab, notch=T, outline=F, ylab="Alternative Allele Frequency (log)", xlab="Stabilising Selection Coefficient")
    boxplot(samp$freq~samp$classSelCoefPos, notch=T, outline=F, ylab="Alternative Allele Frequency (log)", xlab="Directional Selection Coefficient")
    }

  #dev.off()
}


# Site Frequency Spectrum
SNP <- SNP[which(SNP$NbAlt>1),]
neut <- SNP[which(SNP$Conclusion=="Neutral model"),]
neut <- hist(neut$NbAlt/neut$NbTot, breaks=100)
plot(log(neut$counts/sum(neut$counts))~neut$breaks[-1], type="l", xlab="Allele frequency", ylab="log(Density)", main="Site Frequency Spectrum")


stab <- SNP[which(SNP$Conclusion=="Stabilizing model"),]
stab <- hist(stab$NbAlt/stab$NbTot, breaks=100, plot=F)
lines(log(stab$counts/sum(stab$counts))~stab$breaks[-1], col="blue")

pos <- SNP[which(SNP$Conclusion=="Positive model"),]
pos <- hist(pos$NbAlt/pos$NbTot, breaks=100, plot=F)
lines(log(pos$counts/sum(pos$counts))~pos$breaks[-1], col="forestgreen")

legend("topright", legend=c("Neutral", "Stabilising", "Directional"), col=c("black", "blue", "forestgreen"), lty=1, bty="n")
