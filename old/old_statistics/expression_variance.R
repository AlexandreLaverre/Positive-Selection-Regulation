## CTCF analyses
library("ROCR")
library("RColorBrewer") 
library("data.table")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
path <- "/Users/alaverre/Documents/Detecting_positive_selection/"
pathData <- paste0(path, "Tools/JialinTool/")

sp="human"
samp="Wilson"
TFS=c( "CEBPA", "HNF4A", "FOXA1", "HNF6")
minMut = 4
treshold = 0.2
pathResults <- paste(path, "cluster/results/positive_selection/NarrowPeaks", sp, samp, sep="/")

## expression count
geneCount<-fread(paste0(pathData, "data/", sp, "/gtex/getx_liver_gene_tpm.txt"))
colnames(geneCount)[1]<-"geneID"
geneCount$geneID<-gsub("\\..*","",geneCount$geneID)
geneCount[,c(2:176)]<-log2(geneCount[,c(2:176)]+1)
geneCount$var<-apply(geneCount[,c(2:176)],1,function(x) var(x))
geneCount$mean<-apply(geneCount[,c(2:176)],1,function(x) mean(x))


pdf(paste0(path, "results/figures/", sp, "_expression_level_FDR_", treshold, "_", minMut, "mut.pdf"))
for (TF in TFS){
  # Target Genes
  TFBSgene<-fread(paste0(pathResults, "/", TF, "/", TF, ".consensus_peaks.annotatePeaks.txt"))
  colnames(TFBSgene) <- c("peaks_ID", "chr", "start", "end", "TSS_dist", "geneID")
  
  # Plot according to TestPos
  par(mfrow=c(2,2))
  par(mar=c(5,5,4,2))
  
  for (data in c("MLE", "Permutations")){
    if (data == "MLE"){
      TestPos <- read.csv(paste0(pathResults, "/", TF, "/Tests/MLE_summary_exact.csv"))
      TestPos <- TestPos[which(TestPos$Nmut >= minMut),]
      TestPos$FDR <- p.adjust(TestPos$p_value_purif_pos, method="fdr")
      TestPos$Conclusion <- ifelse(TestPos$FDR > treshold , "Neutral model", "Directional model")
      TestPos$pos <- ifelse(grepl("Directional", TestPos$Conclusion), T, F)
      TestPos$null <- ifelse(grepl("Neutral", TestPos$Conclusion), T, F)
    }else{
      TestPos <- fread(paste0(pathResults, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt"))
      TestPos <- TestPos[which(TestPos$NbSub >= minMut),]
      TestPos$FDR <- p.adjust(TestPos$pval.high, method="fdr")
      TestPos$pos <- ifelse(TestPos$FDR<treshold, T, F)
      TestPos$null <- ifelse(TestPos$FDR>=treshold, T, F)
    }
    
    TestPos$peaks_ID <- gsub(".*:\\d+:\\d+[_\\:]", "", TestPos$ID)
    
    Gene <- merge(TestPos, TFBSgene, by="peaks_ID")
    GeneExp <- merge(Gene, geneCount, by="geneID")
    
    ## positive and non-positive sites
    pos <- subset(GeneExp, GeneExp$pos)
    null <- subset(GeneExp, GeneExp$null)
    #posFilter<-pos[!pos$geneID%in%null$geneID,]
    #nullFilter<-null[!null$geneID%in%pos$geneID,]
    
    # Variance
    boxplot(pos$var,null$var,xaxt = "n",
            ylim=c(-0.25,1),notch=T,pch=16,outcex=0.5, main=paste(TF, data), cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
            ylab="Expression variance",col=c(pal[1],pal[2]), outline=F)
    
    text(x=c(1,2),y=-0.36,cex=1.3,srt = 45,adj = 1,labels = c("Positive","Neutral "),xpd = TRUE)
    text(x=c(1,2), y=-0.2, cex=1.3,labels=paste0("n=", c(nrow(pos),nrow(null))))
    
    wtest<-wilcox.test(pos$var,null$var)
    pval<-format(wtest$p.value, scientific = TRUE, digits = 3)
    legend("topleft",legend=paste("p=",pval),bty = 'n',cex=1.5)
    
    # Mean
    boxplot(pos$mean, null$mean, xaxt = "n",
            ylim=c(-0.25,8),notch=T,pch=16,outcex=0.5, main=paste(TF, data), cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
            ylab="Expression mean",col=c(pal[1],pal[2]), outline=F)
    
    text(x=c(1,2),y=-1,cex=1.3,srt = 45,adj = 1,labels = c("Positive","Neutral "),xpd = TRUE)
    
    wtest<-wilcox.test(pos$mean,null$mean)
    pval<-format(wtest$p.value, scientific = TRUE, digits = 3)
    legend("topleft",legend=paste("p=",pval),bty = 'n',cex=1.5)
    
  }
}

dev.off()