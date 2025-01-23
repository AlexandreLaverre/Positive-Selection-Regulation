## CTCF analyses
library("ROCR")
library("RColorBrewer") 
library("data.table")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
path <- "/Users/alaverre/Documents/Detecting_positive_selection/"
pathData <- paste0(path, "Tools/JialinTool/data/")
pathResults <- paste0(path, "results/positive_selection/BroadPeaks/rerun_Jialin_corrected/")
sp="human"

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
#####**** human CTCF adaptive evolution between tissues ****#####
#####*** check the performance of trained model to distinguish binding sites and random sequences ***#####
## five fold cross validation 
predJialin = paste0(pathData, "/", sp, "/CTCF_adaptation/", sp, "_SVM_model/all_merged_ctcf_gkmtrain.cvpred.txt")
predMe = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/CTCF_binding/CTCF/Model/kmer_predicted_weight.txt")
cv<-fread(predJialin)
#cvMe<-fread(predMe)
#predMe <- prediction(cvMe$prediction, cvMe$real_state) 
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
plot(perf,lwd = 3,cex=1.4, main=paste(sp, "CTCF prediction"))
# calculate AUC
auc_result <- performance(pred, measure = "auc")
AUC=unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,paste("AUC =", round(AUC,3)),bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####*** deltaSVM summary ***#####
deltaSVMJialin = paste0(pathData, "/", sp, "/CTCF_adaptation/", sp, "_deltaSVM/ctcf_deltaSVM_highertailTest.txt")
deltaSVMJialin_corrected = paste0(pathResults, "/", sp, "/CTCF_deltaSVM.txt")
deltaSVMMe = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/CTCF_binding/CTCF/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
deltaSVM_ML = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/CTCF_binding/CTCF/Tests/MLE_summary_exact.csv")

deltaSVMMe <- dataMod(fread(deltaSVMMe))
deltaSVMJi <- dataMod(fread(deltaSVMJialin))
deltaSVMJi_cor <- dataMod(fread(deltaSVMJialin_corrected))
deltaSVM_ML <- read.csv(deltaSVM_ML, row.names = 1)

deltaSVM_ML$Conclusion <- ifelse(deltaSVM_ML$p_value_null_pos < 0.01 & deltaSVM_ML$p_value_null_purif > 0.01, "Neutral model", deltaSVM_ML$Conclusion)
deltaSVM_ML$Conclusion <- ifelse(deltaSVM_ML$p_value_null_pos < 0.01 & deltaSVM_ML$p_value_purif_pos < 0.01, "Directional (+)", deltaSVM_ML$Conclusion)

splString<-strsplit(rownames(deltaSVM_ML),"_",fixed=TRUE)
splString<-data.frame(unlist(splString))
ID.bed<-matrix(splString$unlist.splString., ncol=3, byrow=TRUE)[,c(1:3)]
deltaSVM_ML<-cbind(ID.bed, deltaSVM_ML)
colnames(deltaSVM_ML)[c(1:3)]<-c("chr","start","end")
deltaSVM_ML$pValue <- ifelse(grepl("Directional", deltaSVM_ML$Conclusion), 0, 1)

deltaSVMJi_cor <- deltaSVMJi_cor[rownames(deltaSVMJi),]
deltaSVMMe <- deltaSVMMe[rownames(deltaSVMJi),]
deltaSVM_ML <- deltaSVM_ML[rownames(deltaSVMMe),]

list_deltaSVM <- list(deltaSVMJi, deltaSVMJi_cor, deltaSVMMe, deltaSVM_ML)
names(list_deltaSVM) <- c("Jialin_original", "Jialin_corrected", "Mine_Permutation", "Mine_MaxLL")

par(mfrow=c(2,2))
plot(deltaSVMJi$deltaSVM~deltaSVMJi_cor$deltaSVM, cex=0.1, ylim=c(-50, 50), xlim=c(-50, 50), main="DeltaSVM comparison",
     ylab="Permut Jialin", xlab="Permut Jialin corrected")
cor.test(deltaSVMJi$deltaSVM,deltaSVMJi_cor$deltaSVM)

plot(deltaSVMJi$pValue~deltaSVMJi_cor$pValue, cex=0.1, main="pValue comparison",
     ylab="Permut Jialin", xlab="Permut Jialin corrected")
cor.test(deltaSVMJi$pValue, deltaSVMJi_cor$pValue)

plot(deltaSVMJi_cor$deltaSVM~deltaSVMMe$deltaSVM, cex=0.1, ylim=c(-50, 50), xlim=c(-50, 50), main="DeltaSVM comparison",
     ylab="Permut Jialin corrected", xlab="Permut Mine")
cor.test(deltaSVMJi_cor$deltaSVM,deltaSVMMe$deltaSVM)

plot(deltaSVMMe$deltaSVM~deltaSVM_ML$SumObs, cex=0.1, ylim=c(-50, 50), xlim=c(-50, 50), main="DeltaSVM comparison",
     ylab="Permut Mine", xlab="ML (Sum obs)")
cor.test(deltaSVMMe$deltaSVM,deltaSVM_ML$SumObs)

## plot deltaSVM and pvalue
#par(mar=c(7,5,4,2))
n=1.5
par(mfrow=c(1,1))
hist(deltaSVMMe$deltaSVM,breaks = 200,main=" ",xlab="deltaSVM",xlim=c(-30,30),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])

par(mfrow=c(2,2))
for (i in names(list_deltaSVM) ){
  pval=ifelse(i=="Mine_MaxLL", "p_value_null_pos", "pValue")
  hist(list_deltaSVM[[i]][,pval],breaks = 50,xlab="Pvalue",xlim=c(0,1),
       cex.lab=n, cex.axis=n ,col=pal[3], main=i)
}


################################################################################
#####*** positive selection and pleitropy ***#####
organNumb<-fread(paste0(pathData, "/", sp, "/CTCF_adaptation/", sp, "_ctcf_binding/all_merged_annotated.bed"))
## from the fourth column, if the value > 0, indicating this binding site is expressed in this tissue
tissue <- list()
tissue[["human"]] <- c("adrenal gland","B cell","esophagus muscularis mucosa","retinal pigment epithelial cell",
             "omental fat pad","gastrocnemius medialis","astrocyte of the cerebellum","astrocyte of the spinal cord",
             "brain microvascular endothelial cell","choroid plexus epithelial cell","heart left ventricle","liver","upper lobe of left lung",
             "neural cell (in vitro differentiated)","ovary", "pancreas","foreskin fibroblast","peyer patch","prostate gland","lower leg skin",
             "spleen","stomach","testis","thoracic aorta","thyroid gland","tibial nerve","transverse colon","uterus","vagina")

tissue[["mouse"]] <- c("Bone marrow","Cerebellum","Cortical plate","Heart","Kidney","Liver",
                       "Lung","Olfactory bulb","Small intestine","Testis","Thymus")
colnames(organNumb)<-c("chr","start","end", tissue[[sp]])
Ncol=ncol(organNumb)
organNumb$orgNumb<-apply(organNumb[,4:Ncol], 1, function(x) length(x[x>0]))


pdf(file=paste0(path, "results/figures/CTCF_", sp, "_number_organs_new.pdf"))
par(mfrow=c(2,2))
Delta_org <- list()
for (data in names(list_deltaSVM)){
  deltaSVM = list_deltaSVM[[data]]
  deltaSVMorgNumb<-merge(deltaSVM,organNumb,by=c("chr","start","end"))
  posTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue<0.01)
  nonPosTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue>=0.01)
  
  boxplot(posTFBS$orgNumb,nonPosTFBS$orgNumb,ylim=c(0,Ncol),ylab="Number of tissues / cell types",main=data ,xaxt="n",
          col=c(pal[1],pal[2]),notch=T,pch=16,outcex=0.5,cex.lab=1.2,cex.main=1.2)
  text(x=c(1,2),y=-2.5,cex.lab=1.2,srt = 45,adj = 1,labels = c("Positive sites", "Non-positive sites"),xpd = TRUE)
  wTest<-wilcox.test(posTFBS$orgNumb,nonPosTFBS$orgNumb)
  legend("topleft",legend=paste("p=",signif(wTest$p.value,3)),bty = 'n')
  Delta_org <- c(Delta_org, list(deltaSVMorgNumb))
}

names(Delta_org) <- names(list_deltaSVM)
dev.off()

#####*** positive selection and different tissues ***#####
# Define organs to systems and its colors
col <- c("#E41A1C", "#FF7F00", "#66C2A5", "#999999", "#FFFF33", "#F781BF", "#A65628", "#984EA3", "#377EB8", "black", "navy")
names(col) <- c("Nervous", "Male reproductive","Immune", "Endocrine",
                "Integumentary", "Respiratory", "Cardiovascular", "Digestive", 
                "Female reproductive", "Skeletomuscular", "Excretory")
systems <- list()
systems[["human"]] <- c("Endocrine", "Immune", "Digestive", "Nervous", "Integumentary", "Skeletomuscular",
             "Nervous", "Nervous", "Nervous", "Nervous", "Cardiovascular", "Endocrine", 
             "Respiratory", "Nervous", "Female reproductive", "Endocrine", "Male reproductive", "Immune",
             "Male reproductive", "Integumentary", "Immune", "Digestive", "Male reproductive", "Cardiovascular",
             "Endocrine", "Nervous", "Digestive", "Female reproductive", "Female reproductive")

systems[["mouse"]] <- c("Immune", "Nervous", "Nervous", "Cardiovascular", 
                        "Excretory", "Endocrine", "Respiratory", "Nervous", 
                        "Digestive", "Male reproductive", "Immune")

#systems <- paste(systems, "system")
organs <- colnames(organNumb)[c(4:Ncol)]
names(systems[[sp]]) <- organs

pdf(file=paste0(path, "results/figures/CTCF_", sp, "_proportion_organs_new.pdf"), width=8)
for (data in names(list_deltaSVM)){
  deltaSVMorgNumb <- Delta_org[[data]]
  propPos<-c()
  posNumb<-c()
  allNumb<-c()
  for (i in tissue[[sp]]) {
    tissueBasedData<-deltaSVMorgNumb[which(deltaSVMorgNumb[i]>0),]
    propPos[i]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])/nrow(tissueBasedData)
    posNumb[i]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])
    allNumb[i]<-nrow(tissueBasedData)
  }
  Norgan=length(systems[[sp]])
  organProp<-data.frame("organ" = rep(NA,Norgan), "prop" = rep(NA,Norgan),"posNumb" = rep(NA,Norgan),"allNumb" = rep(NA,Norgan))
  organProp$organ<-tissue[[sp]]
  organProp$prop<-propPos
  organProp$posNumb<-posNumb
  organProp$allNumb<-allNumb
  organProp<-organProp[order(organProp$prop),]
  
  ## fisher test for proportions
  ftest<-data.frame("odds_ratio" = rep(NA,Norgan), "pvalue" = rep(NA,Norgan))
  for (i in c(1:Norgan)) {
    temp<-fisher.test(matrix(c(organProp[i,3],organProp[i,4] ,nrow(deltaSVMorgNumb[deltaSVMorgNumb$pValue<0.01,]),nrow(deltaSVMorgNumb)),nrow = 2,ncol = 2))
    ftest[i,1]<-temp$estimate
    ftest[i,2]<-temp$p.value
    
  }
  
  ################################################################################
  ##plot
  par(mfrow=c(1,1))
  par(mar=c(10, 5, 2, 8.5) + 0.1)
  plot(c(1:Norgan), organProp$prop, ylab="Proportion of positive binding sites", main=data,xlab="", 
       pch=16,cex=2,xaxt="n", cex.lab=1.2, cex.main=1.5, cex.axis=1.2,type="p", col=col[systems[[sp]][organProp$organ]],
       ylim=c(min(organProp$prop)-0.005, max(organProp$prop)+0.005))
  
  legend("topright",legend=names(col),col = col, cex=0.9, pch=rep(16,times=10),
         pt.cex=2,bty="n",inset=c(-0.35,0),xpd = TRUE)
  
  abline(v=c(1:Norgan),col="grey", lty=4)
  axis(side = 1, at = c(1:Norgan), labels=F)
  text(c(1:Norgan), par("usr")[3]-0.001, srt = 45,cex=1, adj = 1,labels = organProp$organ, xpd = TRUE)
  ################################################################################
}
dev.off()
