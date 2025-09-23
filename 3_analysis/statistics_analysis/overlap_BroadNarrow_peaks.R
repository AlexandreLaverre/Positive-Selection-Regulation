
path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/peaks_calling/NarrowPeaks/human/Wilson/overlap/"
TFs=c("CEBPA", "FOXA1", "HNF4A", "HNF6")
overlap = list()

for (TF_ref in TFs){
  for (TF_target in TFs){
    if(TF_ref != TF_target){
      ID = paste0(TF_ref, "_to_", TF_target)
      overlap[[ID]] = read.table(paste0(path, ID, ".bed"), h=T)
    }
  }
}

pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/human_peaks_overlap.pdf", width=9)
par(mfrow=c(2,2))
for (ref in TFs){
  mat <- matrix(NA, nrow=3, ncol=2)
  i=1
  evol_ref <- all_MLE_list[[paste0("human_", ref)]]
  for (TF in TFs){
    if(TF != ref){
      ID = paste0(ref, "_to_", TF)
      mat[i,]<- prop.table(table(is.na(overlap[[ID]]$overlap_ID)))
      i=i+1
      
      evol_target <- all_MLE_list[[paste0("human_", TF)]]
      comp = overlap[[ID]]
      rownames(comp) <- comp$ID
      comp <- comp[rownames(evol_ref),]
      comp$ref_evol <- evol_ref$Conclusion
      comp_overlap <- comp[which(!is.na(comp$overlap_ID)),]
      comp_overlap$target_evol <- evol_target[comp_overlap$overlap_ID,]$Conclusion
      
      barplot(prop.table(table(comp_overlap$ref_evol, comp_overlap$target_evol), margin = 2),
              xlab=TF, col=c("forestgreen", "deepskyblue3", "firebrick"), las=1, ylab="Proportion", main=ref)
    }
  }
  rownames(mat) <- TFs[ !TFs == ref]
  colnames(mat) <- c("FALSE", "TRUE")
  barplot(t(mat), col=c("orange", "navy"), main=paste("peaks overlap human", ref))
}
dev.off()


