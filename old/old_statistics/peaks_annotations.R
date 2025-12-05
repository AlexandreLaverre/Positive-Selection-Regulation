library(qvalue)
path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"

species <- c("human","mouse")
TFs <- c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

for (sp in species){
  nb_chr = ifelse(sp =="human", 22, 19)
  chr <- c(seq(1,nb_chr), "X", "Y")
  pdf(paste0(path, "/figures/", sp, "_peaks_annotations.pdf"), width=8)
  for (TF in TFs){
    peaks <- read.table(paste0(path, "positive_selection/", sp, "/new_run/", TF, ".consensus_peaks.annotatePeaks.txt"), h=T, sep="\t")
    peaks_test <- read.table(paste0(path, "positive_selection/", sp, "/new_run/", TF, "_PosSelTest_deltaSVM_10000permutations.txt"), h=T)
    rownames(peaks_test) <- peaks_test$ID
    peaks_signif_ID <- peaks_test[which(qvalue(peaks_test$pval.high)$FDR < 0.01),]$ID
    peaks_test$chr <- sapply(strsplit(peaks_test$ID,":"), `[`, 1)
    peaks_test$start <- sapply(strsplit(peaks_test$ID,":"), `[`, 2)
    peaks_test$end <- sapply(strsplit(peaks_test$ID,":"), `[`, 3)
    
    colnames(peaks)[1] <- "PeakID"
    peaks <- peaks[which(peaks$Chr %in% chr),]
    peaks$Chr <- factor(peaks$Chr, levels=chr)
    peaks$Start <- peaks$Start-1
    peaks$ID <- paste0("chr", peaks$Chr, ":", peaks$Start, ":", peaks$End)
    
    peaks$length <- peaks$End-peaks$Start
    
    # Type
    peaks$type <- sapply(strsplit(peaks$Annotation," "), `[`, 1)
    peaks[which(peaks$type == "promoter-TSS"),]$type <- "promoter"
    peaks$type <- as.factor(peaks$type)
    
    Gene_type <- c("protein_coding", "pseudogene", "lincRNA", "antisense", "miRNA")
    peaks[grep("pseudogene", peaks$Gene.Type),]$Gene.Type <- "pseudogene"
    peaks[which(!peaks$Gene.Type %in% Gene_type),]$Gene.Type <- "other"
    peaks$Gene.Type <- factor(peaks$Gene.Type, levels=c(Gene_type, "other"))
    
    peaks_signif <- peaks[which(peaks$ID %in% peaks_signif_ID),]
    
    par(mfrow=(c(2,2)), mgp=c(2.4,0.6,0))
    par(mai=c(0.6,0.7,0.4,0.2))
    
    # Plot
    hist(peaks$length, xlim=c(0,1500), breaks=100, xlab="Peak length", main=TF)
    hist(peaks$Distance.to.TSS/1000, breaks=300, xlim=c(-150, 150), main="", xlab="Distance to TSS (kb)")
    
    plot(table(peaks$Chr), xlab="Chromosome", ylab="Proportion peaks", las=1)
    plot(table(peaks_signif$Chr)*100/table(peaks$Chr), xlab="Chromosome", ylab="% peaks signif", las=1)
    
    plot(table(peaks$type)*100/nrow(peaks), ylab="% peaks", las=1)
    plot(table(peaks_signif$type)*100/table(peaks$type), ylab="% peaks signif", las=1)


    plot(table(peaks$Gene.Type)*100/nrow(peaks), ylab="% peaks", las=2)
    #plot(table(peaks_signif$Gene.Type)*100/table(peaks$Gene.Type), ylab="Proportion peaks signif", las=1)

    # Gene Lists
    write(peaks[which(peaks$ID %in% peaks_test$ID),]$Entrez.ID, paste0(path, "GOEnrichment/", sp, "/", TF, "_background_genes.txt"))
    write(peaks[which(peaks$ID %in% peaks_signif_ID),]$Entrez.ID, paste0(path, "GOEnrichment/", sp, "/", TF, "_peaks_signif_genes.txt"))
    write.table(peaks_test[,c("chr", "start", "end")], paste0(path, "GOEnrichment/", sp, "/", TF, "_tested_peaks.bed"), row.names=F, col.names=F, sep="\t", quote=F)
    write.table(peaks_test[peaks_signif_ID,c("chr", "start", "end")], paste0(path, "GOEnrichment/", sp, "/", TF, "_signif_peaks.bed"), row.names=F, col.names=F, sep="\t", quote=F)
    
    }
  dev.off()
}
