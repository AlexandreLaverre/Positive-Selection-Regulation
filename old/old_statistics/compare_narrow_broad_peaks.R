path="/Users/alaverre/Documents/Detecting_positive_selection/results/"
TFS=c("HNF6", "HNF4A", "FOXA1", "CEBPA")

pdf(paste0(path, "figures/narrow_vs_broad_peaks_comparisons.pdf"))
for (TF in TFS){
  N_to_B <- read.table(paste0(path, "overlap_narrow_vs_broad/", TF, "_narrow_to_broad.bed"), h=T)
  B_to_N <- read.table(paste0(path, "overlap_narrow_vs_broad/", TF, "_broad_to_narrow.bed"), h=T)
  
  par(mfrow=c(2,2))
  par(mgp=c(2.5, 1.3, 0))
  boxplot(N_to_B$length_frag, B_to_N$length_frag, ylab="Length (bp)", outline=F, notch=T, main=TF, ylim=c(100, 580),
          names=c(paste0("Narrow \n N=", nrow(N_to_B)), paste0("Broad \n N=", nrow(B_to_N))))
  hist(N_to_B$X.overlap, breaks=1, main="Narrow to Broad Peaks", xlab="% overlap (bp)", freq=F)
  hist(B_to_N$X.overlap, breaks=100, main="Broad to Narrow Peaks", xlab="% overlap (bp)", freq=F)
  plot(B_to_N$length_frag, B_to_N$X.overlap, cex=0.1, xlab="Length (bp)", ylab="% overlap", main="Broad to Narrow overlap")
}

dev.off()

