library(data.table)

sp = "mouse"
path <- paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/", sp, "/simulation/")
TFs <- c("CEBPA") #, "FOXA1", "HNF4A")
status <- c("conserved", "specific_loss", "specific_gain")
types <- c("_selection") #""

Test_Pos <- list()
Test_Pos_Sel <- list()
for (TF in TFs){
  for (stat in status){
    for (type in types){
      file = paste0(path, TF, "/calcul_score/PosSelTest_deltaSVM_mouse_triplets_", TF, "_", stat, type, ".txt")
      test <- read.table(file, h=T)
      test$FDR <- p.adjust(test$pval.high, method="fdr")
      
      if (type == ""){Test_Pos[[paste0(TF, "_", stat)]] <- test
      }else{Test_Pos_Sel[[paste0(TF, "_", stat)]] <- test}
      
    }
  }
}

# Delta
DeltaSVM <- lapply(names(Test_Pos), function(x) Test_Pos[[x]]$deltaSVM)
DeltaSVM_simul <- lapply(names(Test_Pos), function(x) Test_Pos[[x]]$med.deltaSVM.simul)
DeltaSVM_Sel_simul <- lapply(names(Test_Pos_Sel), function(x) Test_Pos_Sel[[x]]$med.deltaSVM.simul)

par(mfrow=c(1,3))
boxplot(DeltaSVM, names=names(Test_Pos), notch=T, outline=F, las=1, cex.axis=1.2, cex.lab=1.2,
        xlab="", ylab="delta SVM", main="Observed", col=c("navy", "orange", "forestgreen"))
boxplot(DeltaSVM_simul, names=names(Test_Pos_Sel), notch=T, outline=F, las=1, cex.axis=1.2, cex.lab=1.2,
        xlab="", ylab="delta SVM simul", main="Simul no selection")
boxplot(DeltaSVM_Sel_simul, names=names(Test_Pos_Sel), notch=T, outline=F, las=1, cex.axis=1.2, cex.lab=1.2,
        xlab="", ylab="delta SVM simul", main="Simul with selection")

# Signif
propSignif <- lapply(names(Test_Pos), function(x) length(which(Test_Pos[[x]]$FDR<0.1))/nrow(Test_Pos[[x]]))
propSignif_Sel <- lapply(names(Test_Pos_Sel), function(x) length(which(Test_Pos_Sel[[x]]$FDR<0.1))/nrow(Test_Pos_Sel[[x]]))

par(mfrow=c(1,2))
boxplot(propSignif, names=names(Test_Pos), notch=T, outline=F, las=1, cex.axis=1.2, cex.lab=1.2,
        xlab="", ylab="Prop FDR<0.1", main="Simul no selection", col=c("navy", "orange", "forestgreen"))
boxplot(propSignif_Sel, names=names(Test_Pos_Sel), notch=T, outline=F, las=1, cex.axis=1.2, cex.lab=1.2,
        xlab="", ylab="", main="Simul with selection")
