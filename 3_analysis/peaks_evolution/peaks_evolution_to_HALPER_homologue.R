library(ape)
library(corrplot)
library(RColorBrewer)


# Define paths
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
sp="human"
sample=c("Wilson", "Schmidt12")
TFs=c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

################################################################################
# Retrieve MaxLL results
all_MLE <- list()
for (samp in sample){
  for (TF in TFs){
    MLE.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "_", TF, "_MLE_summary_50bins.csv")
    
    if (file.exists(MLE.file)){
      MLE <- read.csv(MLE.file, h=T, row.names = 1)
      
      #correct error in IDs
      rownames(MLE) <- gsub("\\d+:\\d+:\\d+:", "", rownames(MLE))
      rownames(MLE) <- gsub("_([[:alpha:]])+:\\d+:\\d+:", "_", rownames(MLE))
      MLE$peaks_ID <- gsub(".*:\\d+:\\d+_", "", rownames(MLE))
      
      MLE$Conclusion <- as.factor(MLE$Conclusion)
      MLE$Conclusion <- factor(MLE$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
      MLE$ID <- rownames(MLE)
      
      all_MLE[[TF]] <- MLE
    }
    
    #positive <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Positive model"),]
    #neutral <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Neutral model"),]
    #stabilising <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Stabilizing model"),]
    
    #write(gsub(":", "\t", all_MLE[[TF]]$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_all.bed"))
    #write(gsub(":", "\t", positive$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_positive.bed"))
    #write(gsub(":", "\t", neutral$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_neutral.bed"))
    #write(gsub(":", "\t", stabilising$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_stabilising.bed"))
    }
}

################################################################################
suffix=if (sp == "human") "17way" else "60way.glire"

all_phyloP <- list()
for (TF in TFs){
  phyloP <- read.table(paste0(path, "phyloP/", sp, "/", TF, "_", suffix, ".txt"), h=T)
  rownames(phyloP) <- phyloP$ID
  phyloP <- phyloP[rownames(all_MLE[[TF]]),]
  phyloP$class <- cut(phyloP$Score, quantile(phyloP$Score, probs = seq(0, 1, 0.2),
                                             na.rm=T), include.lowest=T, 
                      labels=c("20%", "40%", "60%", "80%", "100%"))

  all_phyloP[[TF]] <- phyloP
}

################################################################################
par(xpd=T)
pdf(paste0(path, "figures/PhyloP_MaxLLTest_human.pdf"))
for (TF in TFs){
  par(mfrow=c(2,2))
  phyloP <- all_phyloP[[TF]]
  MLE <- all_MLE[[TF]]
  MLE$phyloP <- phyloP$Score
  boxplot(phyloP$Score~MLE$Conclusion, notch=T, outline=F, xlab="Conclusion", ylab="PhyloP score")
  

  prop <- prop.table(table(MLE$Conclusion, phyloP$class), margin = 2)
  barplot(prop, col=c("forestgreen", "deepskyblue2", "firebrick"), las=1, xlab="class phyloP", ylab="Proportion", main=TF)
  legend("top", fill=c("forestgreen", "deepskyblue2", "firebrick"), legend=c("Positive", "Negative", "Neutral"),
         bty="n", ncol=3, y.intersp = -3)
  
  MLE$Conclusion <- NULL
  MLE$ID <- NULL
  CorMat <- cor(MLE, use="pairwise.complete.obs")
  corrplot(CorMat, type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))
}
dev.off()

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

################################################################################
# Retrieve HALPER homologous regions
sp="Homo_sapiens"
targets=c("Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus", "Felis_catus",
          "Canis_lupus_familiaris", "Bos_taurus", "Microcebus_murinus",
          "Manis_javanica", "Pan_troglodytes", "Equus_caballus", 
          "Capra_hircus", "Callithrix_jacchus")

tree <- read.tree("/Users/alaverre/Documents/Detecting_positive_selection/data/species_trees/241-mammals.nk")
# Compute the cophenetic distances
coph_dist <- cophenetic(tree)


pdf(paste0(path, "figures/Distance_from_Homo_sapiens.pdf"), height = 12)
par(mar = c(3, 10, 2, 2))
barplot(sort(coph_dist["Homo_sapiens",]), las=1, main="Genetic Distance from Homo sapiens (ZoonomiaTree)", horiz=T, cex.names=0.6)
dev.off()

DivTree <- c()
for (tar in targets){DivTree <- c(DivTree, coph_dist["Homo_sapiens", tar])}
names(DivTree) <- targets
DivTree <- sort(DivTree)

all_HALPER <- list()
for (TF in TFs){
  HALPER_TF <- list()
  for (target in names(DivTree)){
    HALPER <- read.table(paste0(path, "HALPER/", TF, "/HALPER_", TF, "_", sp, "2", target, ".bed"), h=F)
    colnames(HALPER) <- c("chr", "start", "end", "summit", "peak_ID", "target_length",
                          "peak_length", "summit_to_start", "summit_to_end")
    HALPER[[target]] <- paste(HALPER$chr, HALPER$start, HALPER$end, sep=":")
    HALPER_TF[[target]] <- HALPER[, c("peak_ID", target)]
  }
  
  merged <- Reduce(function(x, y) merge(x, y, by = "peak_ID", all = TRUE), HALPER_TF)
  rownames(merged) <- merged$peak_ID
  all_HALPER[[TF]] <- as.data.frame(merged[,-1])
}

################################################################################
# Plot proportion homologue ~ peaks evolution
par(mfrow=c(2,2))
for (TF in TFs){
  positive <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Positive model"),]
  neutral <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Neutral model"),]
  stabilising <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Stabilizing model"),]
  
  prop_homolog_positive <- apply(all_HALPER[[TF]][positive$peaks_ID,], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_neutral <- apply(all_HALPER[[TF]][neutral$peaks_ID,], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_stabilising <- apply(all_HALPER[[TF]][stabilising$peaks_ID,], 2, function(y) sum(!is.na(y))/length(y))
  
  plot(prop_homolog_neutral~DivTree, ylab="Proportion homologous peak", xlab="Divergence Time", type="b", main=TF, ylim=c(0.4, 1))
  lines(prop_homolog_positive~DivTree, type="b", col="forestgreen")
  lines(prop_homolog_stabilising~DivTree, type="b", col="red")
  legend("topright", legend=c("Stabilising", "Neutral", "Positive"), col=c("red", "black", "forestgreen"), lty=1, bty="n")
}
################################################################################

