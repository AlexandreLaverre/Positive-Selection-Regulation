path = "/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach"
col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Neutral model")

################################################################################
##################### Simulations ##################### 

par(mfrow=c(2,2))
for (distrib in c("Beta", "Gaussian")){
  for (mut in c(10, 20)){
    simul <- read.csv(paste0(path, "/Simulations/MLE_summary_1000simul_100bins_", 
                             mut, "MaxMut_", distrib, ".csv"), h=T)
    simul$Conclusion <- as.factor(simul$Conclusion)
    simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
    simul$Scenario <- as.factor(simul$Scenario)
    simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Stabilizing", "Random"))
    
    # Scenario ~ Conclusion
    table(simul$Scenario, simul$Conclusion)
    
    barplot(table(simul$Conclusion, simul$Scenario), col=col, las=1, cex.lab=1.2,
            xlab="Scenario", ylab="Nb simulation", main=paste(distrib, "MaxMut=", mut))
    legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
    text(x=2, y=1300, labels="Test conclusion", cex=1.2)
    
  }
}


par(mfrow=c(2,2))
for (distrib in c("Beta", "Gaussian")){
  simul <- read.csv(paste0(path, "/Simulations/MLE_summary_1000simul_100bins_", 
                           mut, "MaxMut_", distrib, ".csv"), h=T)
  simul$Conclusion <- as.factor(simul$Conclusion)
  simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
  simul$Scenario <- as.factor(simul$Scenario)
  simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Stabilizing", "Random"))
  
  simul$LRT_null_pos <- (-2 *(simul$LL_neutral - simul$LL_pos))
  simul$p_val_null_pos <- chisq.test(simul$LRT_null_pos)

  for (scenario in c("Random", "Positive", "Stabilizing")){
    hist(simul[which(simul$Scenario==scenario),]$p_value_null_purif, breaks=50, main=paste(distrib, scenario), xlab="Null-Purif")
    hist(simul[which(simul$Scenario==scenario),]$p_value_purif_pos, breaks=50, main=paste(distrib, scenario), xlab="Purif-Pos")
  }

}
    
table(simul$Scenario, simul$Conclusion)



old.gaussian.simul <- read.csv(paste0(path, "/Simulations/MLE_summary_1000simul_100bins.csv"), h=T)
old.gaussian.simul$LRT_null_purif <- -2*(old.gaussian.simul$LL_neutral-old.gaussian.simul$LL_purif)
summary(old.gaussian.simul$LRT_null_purif)
hist(-old.gaussian.simul$LL_neutral-old.gaussian.simul$LL_purif)



pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/ProbaFixEstimation/Comparisons_scenario_conclusion_1000simul_100bins.pdf", width=7, height=8)
par(xpd = TRUE)
par(mfrow=c(2,2))

# Scenario ~ Conclusion
barplot(table(simul$Conclusion, simul$Scenario), col=col, las=1, xlab="Scenario", ylab="Nb simulation", cex.lab=1.2)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
text(x=2, y=1300, labels="Test conclusion", cex=1.2)

# Mean Delta
boxplot(simul$MeanObs~simul$Conclusion*simul$Scenario, notch=T, outline=F,
        names=c("", "Positive", "", "", "Stabilizing", "","", "Random", ""),
        xlab="Scenario", ylab="Mean Delta Obs", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

# Var Delta
boxplot(simul$VarObs~simul$Conclusion*simul$Scenario, notch=F, outline=F,
        names=c("", "Positive", "", "", "Stabilizing", "","", "Random", ""),
        xlab="Scenario", ylab="Var Delta Obs", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

boxplot(simul$Nmut~simul$Conclusion*simul$Scenario, notch=T, outline=F,
        names=c("", "Stabilizing", "", "", "Purification", "","", "Random", ""),
        xlab="Scenario", ylab="Nb subsitution", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
abline(v=c(3.5,6.5), lwd=0.8)

dev.off()

################################################################################
############################ Sequences ######################################### 
species = "human"
TFs = c("CEBPA", "FOXA1", "HNF4A")

obs <- list()
for (TF in TFs){
  data <- read.csv(paste(path, species, TF, "MLE_summary.csv", sep="/"), h=T)
  data$Conclusion <- as.factor(data$Conclusion)
  data$Conclusion <- factor(data$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
  
  pdf(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/human/Stats_", TF, ".pdf"), width=7, height=5)
  par(mfrow=c(1,2))
  table(data$Conclusion)/nrow(data)*100
  barplot(table(data$Conclusion)/nrow(data)*100, col=col, beside=T, main=TF, 
          ylab="% of sequences")
  
  boxplot(log(data$Nmut)~data$Conclusion, main=TF, notch=T, col=col,
          ylab="Substitutions (log)", xlab="Test conclusion", outline=F)
  
  boxplot()
  
  par(mfrow=c(2,2))
  hist(data$pval_NullPurif, breaks=50, main=TF, xlab="pval Null Model vs Purif")
  hist(data$pval_NullPos, breaks=100, main=TF, xlab="pval Null Model vs Positif")
  hist(data$pval_PurifPos, breaks=100, main=TF, xlab="pval Purif vs Positif")
  
  dev.off()
  obs[[TF]] <- data 
  
}

