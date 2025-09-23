library("data.table")
library("gdata")
library("ape")
library("org.Hs.eg.db")
library("ROCR")
library("RColorBrewer") 
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
library("plyr")

#### JIALIN ####
path <- "/Users/alaverre/Documents/Detecting_positive_selection/"
pathData <- paste0(path, "Tools/JialinTool/")
#####*** function: change the format of deltaSVM data ***#####
dataMod<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-1
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:4)])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  return(deltaSVM)
}

deltaSVM<-fread(paste0(pathData, "data/human/deltaSVM/CEBPA/hsap_CEBPA_deltaSVM_highertailTest.txt"))
deltaSVM<-dataMod(deltaSVM)

## target genes of TFBS
TFBSgene<-fread(paste0(pathData, "data/human/TFBS_target_genes/hsap_CEBPA_targetGene.txt"))
TFBSgene<-TFBSgene[,c(1:3,7)]
colnames(TFBSgene)[c(1:4)]<-c("chr","start","end","geneID")
TFBSgene<-unique(TFBSgene)
## expression count
geneCount<-fread(paste0(pathData,"data/human/gtex/getx_liver_gene_tpm.txt"))
colnames(geneCount)[1]<-"geneID"
geneCount$geneID<-gsub("\\..*","",geneCount$geneID)
geneCount[,c(2:176)]<-log2(geneCount[,c(2:176)]+1)
geneCount$var<-apply(geneCount[,c(2:176)],1,function(x) var(x))
geneCount$mean<-apply(geneCount[,c(2:176)],1,function(x) mean(x))
deltaSVMGene<-merge(deltaSVM,TFBSgene,by=c("chr","start","end"))
deltaSVMGeneExp<-merge(deltaSVMGene,geneCount,by="geneID")

## positive and non-positive sites
pos<-subset(deltaSVMGeneExp,deltaSVMGeneExp$pValue<0.01)
nonPos<-subset(deltaSVMGeneExp,deltaSVMGeneExp$pValue>=0.01)
posFilter<-pos[!pos$geneID%in%nonPos$geneID,]
nonPosFilter<-nonPos[!nonPos$geneID%in%pos$geneID,]
pos<-posFilter
nonPos<-nonPosFilter
## plot
par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(pos$var,nonPos$var,xaxt = "n",
        ylim=c(-0.2,2),notch=T,pch=16,outcex=0.5,main="Human CEBPA TFBS",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
        ylab="Expression variance across populations",col=c(pal[1],pal[2]))
text(x=c(1,2),y=-0.36,cex=1.5,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
text(x=c(1,2), y=-0.2, cex=1.5,labels=paste0("n=", c(nrow(pos),nrow(nonPos))))

wtest<-wilcox.test(pos$var,nonPos$var)
legend("topleft",legend=paste("p=",signif(wtest$p.value,3)),bty = 'n',cex=1.5)

#### ME ####
sp= "human"
TF="CEBPA"

peaks_test <- read.table(paste0(path, "results/positive_selection/", sp, "/new_run/", TF, "_PosSelTest_deltaSVM_10000permutations.txt"), h=T)
rownames(peaks_test) <- peaks_test$ID
peaks_test$FDR <- p.adjust(peaks_test$pval.high, method="fdr")

peaks <- read.table(paste0(path, "results/positive_selection/", sp, "/new_run/", TF, ".consensus_peaks.annotatePeaks.txt"), h=T, sep="\t")
colnames(peaks)[1] <- "PeakID"
peaks <- peaks[which(peaks$Chr %in% chr),]
peaks$Chr <- factor(peaks$Chr, levels=chr)
peaks$Start <- peaks$Start-1
peaks$ID <- paste0("chr", peaks$Chr, ":", peaks$Start, ":", peaks$End)
peaks2gene <- as.data.frame(cbind(peaks$ID, peaks$Entrez.ID))
colnames(peaks2gene) <- c("ID", "geneID")

peaks2gene <- peaks2gene[which(peaks2gene$ID %in% peaks_test$ID),]
deltaSVMGene_me<-merge(peaks_test,peaks2gene,by=c("ID"))
deltaSVMGeneExp_me<-merge(deltaSVMGene_me,geneCount,by="geneID")

pos<-subset(deltaSVMGeneExp_me,deltaSVMGeneExp_me$FDR<=0.1)
nonPos<-subset(deltaSVMGeneExp_me,deltaSVMGeneExp_me$FDR>0.1)
posFilter<-pos[!pos$geneID%in%nonPos$geneID,]
nonPosFilter<-nonPos[!nonPos$geneID%in%pos$geneID,]
pos<-posFilter
nonPos<-nonPosFilter

