library(ROCR)

sp="human"
TF="Wilson/HNF4A"
peakType="NarrowPeaks"
path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/", peakType, "/", sp, "/", TF, "/")
pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

bins = c("exact_50", "exact_absolute_50")
datas = c("beta_0.0Null_25Stab_25Pos", "beta_0.1Null_25Stab_25Pos", "beta_0.25Null_25Stab_25Pos")
selection = c("positive", "neutral")

col=c("#4393C3", "#D6604D")
names(col)=bins

#### Find treshold point in predictions #### 
find_threshold_point <- function(pred, threshold) {
  fpr <- unlist(slot(performance(pred, "fpr"), "y.values"))
  tpr <- unlist(slot(performance(pred, "tpr"), "y.values"))
  scores <- unlist(slot(pred, "predictions"))
  
  # Find the index of the closest threshold
  idx <- which.min(abs(scores - threshold))
  
  return(c(fpr[idx], tpr[idx]))
}

#### Check each simul and test ####
for (data in datas){
  for (bin in bins){
    simul_type <- list()
    for (sel in selection){
      simul <- read.csv(paste0(path, "/Tests/MLE_summary_simulated_", data, "_", sel, "_", bin, "bins_threshold_0.01.csv"), h=T, row.names = 1)
      permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_", data, "_", sel, ".txt"), h=T, row.names = 1)
      permut <- permut[rownames(simul),]
      
      simul <- cbind(simul, permut[,c("pval.high", "deltaSVM")])
      simul$real_state <- ifelse(sel == "positive", 1, 0)
      
      # Convert p-values to a confidence score (higher means more confidence)
      simul$score <- -log10(simul$p_value_null_pos) 
      simul$pval.permut <- ifelse(simul$deltaSVM>0, simul$pval.high, 1-simul$pval.high)
      simul$score_permut <- -log10(simul$pval.permut + runif(nrow(simul), 1e-7, 1e-6)) #add small noise to avoid 0
      
      simul_type[[sel]] <- simul
    }
    
    all_simul <- rbind(simul_type[["positive"]], simul_type[["neutral"]])

    # Initiate ROC plot with Permutations Test
    if (bin == "exact_50"){
      # Get prediction and performance
      pred_permut <- prediction(all_simul$score_permut, all_simul$real_state)
      perf_permut <- performance(pred_permut, "tpr", "fpr" )
      pval.permut.1 <- find_threshold_point(pred_permut, -log10(0.01))
      pval.permut.5 <- find_threshold_point(pred_permut, -log10(0.05))
      
      plot(perf_permut, lwd = 2,cex=1.2, main=data)
      points(pval.permut.1[1], pval.permut.1[2], pch = 19, cex = 1.5) 
      points(pval.permut.5[1], pval.permut.5[2], pch = 1, cex = 1.5) 
      
      # add AUC value
      auc <- performance(pred_permut, measure = "auc")
      AUC=unlist(slot(auc, "y.values"))
      legend(0.5,0.1, paste("AUC permut. =", round(AUC,3)),bty="n",lwd = 3, cex=1.1) 
      
    }
    
    # Get prediction and performance
    pred_MLL <- prediction(all_simul$score, all_simul$real_state)
    perf_MLL <- performance(pred_MLL, "tpr", "fpr" )
    pval.MLL.1 <- find_threshold_point(pred_MLL, -log10(0.01)) 
    pval.MLL.5 <- find_threshold_point(pred_MLL, -log10(0.05)) 
    
    plot(perf_MLL, col=col[bin], lwd = 2, cex=1.2, add=T)
    points(pval.MLL.1[1], pval.MLL.1[2], col = col[bin], pch = 19, cex = 1.5)  # Ranked test
    points(pval.MLL.5[1], pval.MLL.5[2], col = col[bin], pch = 1, cex = 1.5)  # Ranked test
    
    auc <- performance(pred_MLL, measure = "auc")
    AUC = unlist(slot(auc, "y.values"))
    if (bin == "exact_50"){
      legend(0.5, 0.3, paste("AUC MLL rank.=", round(AUC,3)),bty="n",lwd = 3, cex=1.1, col = col[bin]) 
    }else{
      legend(0.5, 0.2, paste("AUC MLL abs.=", round(AUC,3)),bty="n",lwd = 3, cex=1.1, col = col[bin]) 
    }

  }
}
