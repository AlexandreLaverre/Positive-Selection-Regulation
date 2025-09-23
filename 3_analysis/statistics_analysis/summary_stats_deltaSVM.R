###########################################################################
#####**** Libraries and functions ****#####
library("data.table")
library("gdata")
library("ape")
library("org.Hs.eg.db")
library("ROCR")
library("RColorBrewer")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
library("plyr")
library("qvalue")

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

###########################################################################
#####**** Paths ****#####
args = commandArgs(trailingOnly=TRUE)
sp = args[1] #"human"
sample = args[2] #"CEBPA"
from = args[3] # "Jialin" or "Mine"

focal = ifelse(sp=="mouse", "C57BL/6J", sp)

path = "/Users/alaverre/Documents/Detecting_positive_selection/"
pathResults = ifelse(from == "Jialin", paste0(path, "Tools/JialinTool/data/", sp),
                     paste0(path,"results/positive_selection/rerun_Jialin_corrected/", sp))

###########################################################################
#####**** Data ****#####
# SVM model
model = list.files(path=paste0(path, "Tools/JialinTool/data/", sp, "/SVM_model/", sample, "/"), pattern="*cvpred.txt", full.names=T)
cv <- fread(model)
colnames(cv)<-c("position","prediction","real_state","cv_number")

# Delta SVM
if (sp=="mouse"){states = c("conserved", "gain", "loss")
}else{states = c("_")}

deltaSVM <- list()
for (state in states){
  file = list.files(path=paste0(pathResults, "/deltaSVM/", sample, "/"), pattern=state, full.names=T)[1]
  deltaSVM[[state]] <- dataMod(fread(file))
  if (state == "loss" & from == "Mine"){deltaSVM[[state]]$pValue <- 1-deltaSVM[[state]]$pValue}
  
  ## pvalue correction for multiple tests
  deltaSVM[[state]]$FDR<-p.adjust(deltaSVM[[state]]$pValue, method="fdr")
  #deltaSVM[[state]]$FDR <- qvalue(deltaSVM[[state]]$pValue)$qvalues
}

# Binding
if (sp =="mouse"){
  ChIP = list.files(path=paste0(path, "Tools/JialinTool/data/", sp, "/ChIP-Seq/"), pattern=paste0(sample,".txt"), full.names=T)
  binding_intensity<-fread(ChIP)
  binding_intensity$chr<-paste0(rep("chr",nrow(binding_intensity)),binding_intensity$chr)
  colnames(binding_intensity)[c(2,3)]<-c("start","end")
}


###########################################################################
#####**** Figures ****#####
pdf(paste0(path, "results/figures/", sp, "_", sample, "_", from, "_qvalue_deltaSVM_summary.pdf"), width=8)
##** 1 - Model performance **##
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )

plot(perf,lwd = 3,cex=1.4, main=paste(sp, sample, "prediction"))

auc_result <-performance( pred, measure = "auc")
AUC <- signif(unlist(slot(auc_result, "y.values")), 3)
legend(0.3,0.6, paste0("(AUC = ",  AUC, ")"), bty="n",lwd = 3,cex=1.4) 

##** 2 - deltaSVM summary **##
par(mfrow=c(3,3))
par(mar=c(5,5,2,1))
CEX=1.5

for (state in states){
  delta <- deltaSVM[[state]]
  
  hist(delta$deltaSVM, breaks = 60, main="", xlab="deltaSVM", xlim=c(-30,30),
       cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[4])
  hist(delta$pValue, breaks = 60, main=state, xlab="p-value", xlim=c(0,1),
       cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[3])
  hist(delta$FDR, breaks = 60, main="", xlab="q-value", xlim=c(0,1),
       cex.lab=CEX, cex.axis=CEX, cex.main=2, col=pal[3])
}

##** 3 - Proportion of positive selection **##
par(mfrow=c(1,2))
par(mar=c(8,6,8,0.5))
PosProp_FDR <- c()
PosProp_pval <- c()

for (state in states){
  delta <- deltaSVM[[state]]
  Pos = nrow(delta[which(delta$pValue<0.01),])
  PosProp_pval[state] <- Pos/nrow(delta)
  
  Pos_FDR = nrow(delta[which(delta$FDR<0.1),])
  PosProp_FDR[state] <- Pos_FDR/nrow(delta)
}

barplot(PosProp_pval, col=pal[2], cex.lab=CEX, cex.main=CEX, cex.axis=CEX, cex.names = CEX, las=3,
        ylab="Proportion of positive binding sites \n (pval<0.01)", main=paste(sample, "binding sites"))

barplot(PosProp_FDR, col=pal[2], cex.lab=CEX, cex.main=CEX, cex.axis=CEX, cex.names = CEX, las=3,
        ylab="(FDR<0.1)", main="")

##** 4 - Delta SVM and binding intensity **##
if (sp == "mouse"){
  par(mfrow=c(1,2))
  par(mar=c(9,6,4,1))
  CEX=1.2
  for (test in c("pValue", "FDR")){
    treshold = ifelse(test=="Pvalue", 0.01, 0.1)
    
    for (state in c("conserved", "gain")){ # cannot be evaluated in loss
      deltaSVM_intensity <- merge(deltaSVM[[state]], binding_intensity, by=c("chr","start","end"))
      pos <-subset(deltaSVM_intensity, deltaSVM_intensity[test]<treshold)
      nonPos<-subset(deltaSVM_intensity, deltaSVM_intensity[test]>=treshold)
      
      ylab=ifelse(state=="conserved", paste("Binding intensity in", focal, "\n", test, "<", treshold),"") #only in first plot
      boxplot(pos$intensity, nonPos$intensity, col=c(pal[1],pal[2]),
              ylim= c(-50,1500), ylab=ylab, cex.axis=CEX,
              notch=T, pch=16, outcex=0.5, main=state, cex.lab=CEX, cex.main=2)
      
      text(x=c(1,2), y=-200, cex=CEX, srt = 45, adj = 1, xpd = TRUE,
           labels = c("Positive sites","Non-positive sites"))
      
      tryCatch({
        wtest <- wilcox.test(pos$intensity,nonPos$intensity)
        pval <- signif(wtest$p.value,3)
      }, error=function(e){pval <<- "NA"})
      
      text(x=1.5, y=1600, labels=paste("pval=",pval), cex=CEX, xpd = TRUE,)
      text(x=c(1,2), y=-50, labels=paste0("n=", c(nrow(pos),nrow(nonPos))), cex=CEX)
    }
  }
  
}

dev.off()
###########################################################################
