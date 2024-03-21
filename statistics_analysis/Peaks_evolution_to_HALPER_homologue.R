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
    MLE.file <- paste0(path, "MaxLikelihoodApproach/real_peaks/", sp, "_", samp, "_", TF, "_MLE_summary_50bins.csv")
    
    if (file.exists(MLE.file)){
      MLE <- read.csv(MLE.file, h=T, row.names = 1)
      MLE$Conclusion <- as.factor(MLE$Conclusion)
      MLE$Conclusion <- factor(MLE$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
      MLE$ID <- rownames(MLE)
      MLE$species <- sp
      
      all_MLE[[TF]] <- MLE
    }}}

################################################################################
# Retrieve HALPER homologous regions
sp="Homo_sapiens"
targets=c("Macaca_mulatta", "Mus_musculus", "Mus_spretus", "Mus_caroli", "Felis_catus",
          "Canis_lupus_familiaris", "Oryctolagus_cuniculus", "Rattus_norvegicus", "Bos_taurus")

all_HALPER <- list()
for (TF in TFs){
  HALPER_TF <- list()
  for (target in targets){
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
  
  positive <- rownames(all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Positive model"),])
  neutral <- rownames(all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Neutral model"),])
  stabilising <- rownames(all_MLE[[TF]][which(all_MLE[[TF]]$Conclusion == "Stabilizing model"),])
  
  prop_homolog_positive <- apply(all_HALPER[[TF]][positive,], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_neutral <- apply(all_HALPER[[TF]][neutral,], 2, function(y) sum(!is.na(y))/length(y))
  prop_homolog_stabilising <- apply(all_HALPER[[TF]][stabilising,], 2, function(y) sum(!is.na(y))/length(y))
  
  plot(prop_homolog_neutral, ylab="Proportion homologous peak", xlab="Divergence Time", type="b", main=TF, ylim=c(0.5, 1))
  lines(prop_homolog_positive, type="b", col="forestgreen")
  lines(prop_homolog_stabilising, type="b", col="red")
  legend("topright", legend=c("Negative", "Neutral", "Positive"), col=c("red", "black", "forestgreen"), lty=1, bty="n")
}
################################################################################

