path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/human/"
TFS= c("Wilson/CEBPA", "Wilson/FOXA1", "Wilson/HNF4A", "Wilson/HNF6", "Schmidt12/CTCF") 

for (TF in TFS){
  message(TF)
  peaks_overlap <- read.table(paste0(path, TF, "/peaks_overlap_transposons.bed"), h=T)
  rownames(peaks_overlap) <- peaks_overlap$end
  
  peaks_evol <- read.csv(paste0(path, TF, "/MLE_summary_quantile_50bins_threshold_0.01.csv"), h=T)
  peaks_evol$short_ID <- sub("_.*", "", peaks_evol$ID)
  peaks_evol <- peaks_evol[grepl(":", peaks_evol$short_ID, ignore.case = TRUE), ]
  
  rownames(peaks_evol) <- peaks_evol$short_ID
  peaks_evol$Transposon <- peaks_overlap[rownames(peaks_evol),]$X.overlap
  peaks_evol$PA_Transposon <- ifelse(peaks_evol$Transposon>0, "with TE", "w.o TE")
  
  boxplot(peaks_evol$Transposon~peaks_evol$Conclusion, xlab="", ylab="% Transposon", main=TF)
  barplot(prop.table(table(peaks_evol$PA_Transposon,peaks_evol$Conclusion), margin = 2), 
          ylab="Proportion", legend.text = TRUE, col = c("skyblue", "salmon"), main=TF)
  abline(h=0.5, cex=3)
}
