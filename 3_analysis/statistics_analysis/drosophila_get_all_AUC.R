library("ROCR")
library("data.table")
library("progress")

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
exp.folders <- list.dirs(path, recursive = FALSE)

# All results
method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, ".Rds")
all_MLE_list <- readRDS(AllMerged)

# Get Model performances
all_AUC <- data.frame(exp = rep(NA,length(exp.folders)), AUC = NA, SampleSize=NA)
pb <- progress_bar$new(format = "  Processing [:bar] :percent (:elapsed s) | ETA :eta", total = length(exp.folders), clear = FALSE, width = 60)
for (i in 1:length(exp.folders)){
  pb$tick()
  exp.path = exp.folders[i]
  exp <- basename(exp.path)
  if (exp == "log"){next}
  
  cv<-fread(paste0(exp.path, "/Model/", exp, ".cvpred.txt"))
  colnames(cv)<-c("position","prediction","real_state","cv_number")
  pred <- prediction(cv$prediction, cv$real_state) 
  auc_result <- performance(pred, measure = "auc")
  AUC = signif(unlist(slot(auc_result, "y.values")), digits=3)
  
  all_AUC$exp[i] <- exp
  all_AUC$AUC[i] <- AUC
  all_AUC$SampleSize[i] <- unique(all_MLE_list[[exp]]$SampleSize)
}

all_AUC <- all_AUC[!is.na(all_AUC$AUC),]
rownames(all_AUC) <- all_AUC$exp

# Plots
hist(all_AUC$AUC, breaks=100, xlab="AUC", main="Drosophila - all datasets")
plot(all_AUC$AUC~all_AUC$SampleSize, cex=0.8, xlab="Sample Size", ylab="AUC", pch=19)
abline(v=1500, col="red")



