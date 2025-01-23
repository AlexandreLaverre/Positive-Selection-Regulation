##### Polymorphism vs Substitutions
library("ROCR")
library("RColorBrewer") 
library("data.table")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
path <- "/Users/alaverre/Documents/Detecting_positive_selection/"
pathData <- paste0(path, "Tools/JialinTool/data/")
pathResults <- paste0(path, "results/positive_selection/BroadPeaks/rerun_Jialin_corrected/")

dataMod<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  if (ncol(deltaSVM) >= 7){
    deltaSVM <- deltaSVM[,c("ID", "deltaSVM", "NbSub", "pval.high")]
    nsplit=3;minus=0
  }else{nsplit=5; minus=1}
  
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=nsplit, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-minus
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:ncol(deltaSVM))])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$ID <- paste(ID.bed[,1],ID.bed[,2],ID.bed[,3], sep="_")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  
  rownames(deltaSVM) <- deltaSVM$ID
  return(deltaSVM)
}

################################################################################
TFS = c("CEBPA", "HNF4A")
sp = "human"
pdf(file=paste0(path, "results/figures/", sp, "_substitutions_polymorphism.pdf"), width=8)
par(mfrow=c(2,2))
par(mar=c(9,5,4,2))

for (TF in TFS){
  deltaSVMJialin = paste0(pathData, "/", sp, "/deltaSVM/", TF, "/deltaSVM_highertailTest.txt")
  deltaSVMJialin_corrected = paste0(pathResults, "/", sp, "/deltaSVM/", TF, "/", TF, "_deltaSVM.txt")
  deltaSVMMe = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/Wilson/",TF , "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
  deltaSVM_ML = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/Wilson/", TF, "/Tests/MLE_summary_exact.csv")
  
  deltaSVMJi <- dataMod(fread(deltaSVMJialin))
  deltaSVMJi_cor <- dataMod(fread(deltaSVMJialin_corrected))
  deltaSVMMe <- fread(deltaSVMMe)
  deltaSVM_ML <- read.csv(deltaSVM_ML, row.names = 1)
  deltaSVM_ML$pValue <- ifelse(grepl("Directional", deltaSVM_ML$Conclusion), 0, 1)
  
  # Poly
  HumanPolymorphism <- fread(paste0(path, "/cluster/results/polymorphism_analyses/NarrowPeaks/", sp, "/Wilson/", TF, "/VCF/overlap_peaks.vcf"))
  NbPoly <- as.data.frame(table(HumanPolymorphism$V2))
  rownames(NbPoly) <- NbPoly$Var1
  deltaSVMMe$NbPoly <- NbPoly[deltaSVMMe$ID,]$Freq
  deltaSVM_ML$NbPoly <- NbPoly[rownames(deltaSVM_ML),]$Freq
  
  JialinPolymorphism <- fread(paste0(pathData, "/", sp, "/substitutions_polymorphisms/", TF, "_sub_poly.txt"))
  deltaSVMJi<-merge(deltaSVMJi,JialinPolymorphism,by=c("chr","start","end"))
  deltaSVMJi_cor<-merge(deltaSVMJi_cor,JialinPolymorphism,by=c("chr","start","end"))
  
  #homogeneous
  deltaSVM_ML$NbSub <- deltaSVM_ML$Nmut
  deltaSVMJi$NbSub <- deltaSVMJi$subNumb
  deltaSVMJi$NbPoly <- as.integer(deltaSVMJi$polyNumb)
  deltaSVMJi_cor$NbSub <- deltaSVMJi_cor$subNumb
  deltaSVMJi_cor$NbPoly <- as.integer(deltaSVMJi_cor$polyNumb)
  deltaSVMMe$pValue <- deltaSVMMe$pval.high
  
  # cor
  cor.test(deltaSVMMe$NbSub,as.integer(deltaSVMMe$NbPoly))
  cor.test(deltaSVM_ML$NbSub,as.integer(deltaSVM_ML$NbPoly))
  cor.test(deltaSVMJi$NbSub,as.integer(deltaSVMJi$NbPoly))
  cor.test(deltaSVMJi_cor$NbSub,as.integer(deltaSVMJi_cor$NbPoly))
  
  ## positive and non-positive sites
  list_deltaSVM <- list(deltaSVMJi, deltaSVMJi_cor, deltaSVMMe, deltaSVM_ML)
  names(list_deltaSVM) <- c("Jialin_original", "Jialin_corrected", "Mine_Permutation", "Mine_MaxLL")
  
  for (data in names(list_deltaSVM)){
    deltaSVM = as.data.frame(na.omit(list_deltaSVM[[data]]))
    
    pos<- deltaSVM[which(deltaSVM$pValue<0.01),]
    nonPos<- deltaSVM[which(deltaSVM$pValue>=0.01),]
    subNumb<-c(sum(pos$NbSub),sum(nonPos$NbSub))
    polyNumb<-c(sum(pos$NbPoly),sum(nonPos$NbPoly))
    
    if (data %in% c("Jialin_original", "Jialin_corrected")){ymax=1.5}else{ymax=0.35}
    ## plot
    bp<-barplot(subNumb/polyNumb,ylim=c(0,ymax),main=paste0(data, " (", TF, ")"), cex.lab=1.2,cex.main=1.5,cex.axis =1.2,ylab="# Substitutions / # polymorphisms",col=c(pal[1],pal[2]))
    text(x=bp,y=0-1.5/15,cex=1.5,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
    text(x=bp, y=0+1.5/15, labels=paste0("n=", c(nrow(pos),nrow(nonPos))),cex=1.5)
    # fisher exact test
    fTest<-fisher.test(matrix(c(subNumb,polyNumb),nrow = 2,ncol = 2))
    legend("topleft",legend=paste("p=",signif(fTest$p.value,3)),bty = 'n',cex=1.5)
  }
}

dev.off()
