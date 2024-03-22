library(ape)

# Define paths
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
sp="mouse"
sample=c("Wilson", "Schmidt12")
TFs=c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

################################################################################
# Retrieve MaxLL results
all_MLE <- list()
for (samp in sample){
  for (TF in TFs){
    MLE.file <- paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "_", samp, "_", TF, "_MLE_summary_50bins.csv")
    
    if (file.exists(MLE.file)){
      MLE <- read.csv(MLE.file, h=T, row.names = 1)
      MLE$Conclusion <- as.factor(MLE$Conclusion)
      MLE$Conclusion <- factor(MLE$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
      MLE$ID <- rownames(MLE)
      MLE$species <- sp
      
      all_MLE[[TF]] <- MLE
    }
    
    positive <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Positive model"),]
    neutral <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Neutral model"),]
    stabilising <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Stabilizing model"),]
    
    write(gsub(":", "\t", all_MLE[[TF]]$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_all.bed"))
    write(gsub(":", "\t", positive$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_positive.bed"))
    write(gsub(":", "\t", neutral$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_neutral.bed"))
    write(gsub(":", "\t", stabilising$ID), file=paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "/", TF, "_stabilising.bed"))
    }
}


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
  ###################### Attribute peaks ID: TEMPORARY!!! ######################
  all_MLE[[TF]]$peak_ID <- sample(rownames(all_HALPER[[TF]]), nrow(all_MLE[[TF]]), replace=F)
  rownames(all_MLE[[TF]]) <- all_MLE[[TF]]$peak_ID
  ##############################################################################
  
  positive <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Positive model"),]
  neutral <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Neutral model"),]
  stabilising <- all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Stabilizing model"),]
  
  prop_homolog_positive <- apply(all_HALPER[[TF]][rownames(positive),], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_neutral <- apply(all_HALPER[[TF]][rownames(neutral),], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_stabilising <- apply(all_HALPER[[TF]][rownames(stabilising),], 2, function(y) sum(!is.na(y))/length(y))
  
  plot(prop_homolog_neutral~DivTree, ylab="Proportion homologous peak", xlab="Divergence Time", type="b", main=TF, ylim=c(0.5, 1))
  lines(prop_homolog_positive~DivTree, type="b", col="forestgreen")
  lines(prop_homolog_stabilising~DivTree, type="b", col="red")
  legend("topright", legend=c("Negative", "Neutral", "Positive"), col=c("red", "black", "forestgreen"), lty=1, bty="n")
}
################################################################################

