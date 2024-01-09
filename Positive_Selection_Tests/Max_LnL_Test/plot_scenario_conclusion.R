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

################################################################################
################################################################################
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

################################################################################
################################################################################
# Simulations 500 sequences
library(corrplot)
library(RColorBrewer)

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Neutral model")

bin = "50"
datas = c("", "_500simul")
selection = c("neutral", "stabilising", "positive")
maxSub=150
par(mfrow=c(2,2))

all_delta <- list()
all_simul <- list()

for (data in datas){
  for (sel in selection){
  obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
  deltas <- read.table(paste0(path, "/positive_selection/human/Wilson/CEBPA/deltas/simulated", data, "_", sel, "_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
  row.names(deltas) <- deltas$seq_name
  all_delta[[paste0(sel, data)]] <- deltas
  
  simul <- read.csv(paste0(path, "MaxLikelihoodApproach/human/CEBPA/simulations/MLE_summary_simulated_", sel, "_", bin, "bins", data, ".csv"), h=T, row.names = 1)
  simul$Conclusion <- as.factor(simul$Conclusion)
  simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
  
  simul$deltaSVM <- deltas[row.names(simul),]$deltaSVM
  simul$Diffdelta <-simul$deltaSVM-simul$SumObs
  hist(simul$Diffdelta, breaks=20, xlab="diff deltaSVM", main=paste(sel, data))
  all_simul[[paste0(sel, data)]] <- simul
  
  print(paste(sel, data))
  print(table(simul$Conclusion))
  }
  
  
}

simul$Conclusion <- NULL
CorMat <- cor(simul)
corrplot(CorMat, type="upper", order="hclust", main=type, col=brewer.pal(n=8, name="RdYlBu"))

for (data in datas){
  rand <- all_simul[[paste0("neutral", data)]]
  rand$ID <- rownames(rand)
  rand$Scenario <- "Random"
  rownames(rand) <- NULL
  pos <- all_simul[[paste0("positive", data)]]
  pos$ID <- rownames(pos)
  pos$Scenario <- "Positive"
  rownames(pos) <- NULL
  stab <- all_simul[[paste0("stabilising", data)]]
  stab$ID <- rownames(stab)
  stab$Scenario <- "Stabilizing"
  rownames(stab) <- NULL
  simul <- rbind(rand, pos, stab)
  simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Stabilizing", "Random"))
  
  boxplot(simul$Nmut~simul$Scenario, notch=T)
  boxplot(simul$Nmut~simul$Conclusion, notch=T)
  boxplot(simul$Diffdelta~simul$Scenario, notch=T, outline=F)
  boxplot(simul$Diffdelta~simul$Conclusion, notch=T, outline=F)
  
  
  pdf(paste0(path, "/figures/Simulations", data, "_epistatic_effect.pdf"))
  par(mfrow=c(2,2))
  # Barplot proportion epistasis per substitution number
  simul$Epistasis <- ifelse(simul$Diffdelta==0, 0, 1)
  barplot(tapply(simul$Epistasis, simul$Nmut, function(x) (sum(x)/length(x))*100),
          xlab="Number of substitution", ylab="% seq with epistasis", main="All simulations", las=1)
  
  barplot(tapply(simul[which(simul$Scenario=="Random"),]$Epistasis, simul[which(simul$Scenario=="Random"),]$Nmut, function(x) (sum(x)/length(x))*100),
          xlab="Number of substitution", ylab="% seq with epistasis", main="Random simulations", las=1)
  
  barplot(tapply(simul[which(simul$Scenario=="Positive"),]$Epistasis, simul[which(simul$Scenario=="Positive"),]$Nmut, function(x) (sum(x)/length(x))*100),
          xlab="Number of substitution", ylab="% seq with epistasis", main="Positive simulations", las=1)
  
  barplot(tapply(simul[which(simul$Scenario=="Stabilizing"),]$Epistasis, simul[which(simul$Scenario=="Stabilizing"),]$Nmut, function(x) (sum(x)/length(x))*100),
          xlab="Number of substitution", ylab="% seq with epistasis", main="Stabilising simulations", las=1)
  
  # Boxplot epistatic effect
  boxplot(simul$Diffdelta~simul$Nmut, notch=T, outline=F, xlab="Number of substitution", ylab="Epistatic effect", main="All simulations")
  abline(h=0, lty="dotted", cex=1.2)
  prop=1-length(which(simul$Diffdelta == 0))/nrow(simul)
  mtext(paste("Proportion of epistasis =", signif(prop, digits=2)))
  
  boxplot(rand$Diffdelta~rand$Nmut, notch=T, outline=F, xlab="Number of substitution", ylab="Epistatic effect", main="Random simulations")
  abline(h=0, lty="dotted", cex=1.2)
  prop=1-length(which(rand$Diffdelta == 0))/nrow(rand)
  mtext(paste("Proportion of epistasis =", signif(prop, digits=2)))
  
  boxplot(pos$Diffdelta~pos$Nmut, notch=T, outline=F, xlab="Number of substitution", ylab="Epistatic effect", main="Positive simulations")
  abline(h=0, lty="dotted", cex=1.2)
  prop=1-length(which(pos$Diffdelta == 0))/nrow(pos)
  mtext(paste("Proportion of epistasis =", signif(prop, digits=2)))
  
  boxplot(stab$Diffdelta~stab$Nmut, notch=T, outline=F, xlab="Number of substitution", ylab="Epistatic effect", main="Stabilising simulations")
  abline(h=0, lty="dotted", cex=1.2)
  prop=1-length(which(stab$Diffdelta == 0))/nrow(stab)
  mtext(paste("Proportion of epistasis =", signif(prop, digits=2)))
  
  dev.off()
}


barplot(table(simul$Conclusion, simul$Scenario), col=col, las=1, cex.lab=1.2, xlab="Scenario", ylab="Nb simulation")
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
text(x=2, y=1300, labels="Test conclusion", cex=1.2)