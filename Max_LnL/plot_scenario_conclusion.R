path = "/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach"
col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Neutral model")

################################################################################
##################### Simulations ##################### 
simul <- read.csv(paste0(path, "/ProbaFixEstimation/MLE_summary_100simul_30bins.csv"), h=T)
simul$Conclusion <- as.factor(simul$Conclusion)
simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
simul$Scenario <- as.factor(simul$Scenario)
simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Stabilizing", "Random"))

table(simul$Scenario, simul$Conclusion)

pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/ProbaFixEstimation/Comparisons_scenario_conclusion.pdf", width=7, height=8)
par(xpd = TRUE)
par(mfrow=c(2,1))

# Scenario ~ Conclusion
barplot(table(simul$Conclusion, simul$Scenario), col=col, las=1, xlab="Scenario", ylab="Nb simulation", cex.lab=1.2)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
text(x=2, y=1300, labels="Test conclusion", cex=1.2)

# Mean Delta
boxplot(simul$MeanObs~simul$Conclusion*simul$Scenario, notch=T, outline=F,
        names=c("", "Positive", "", "", "Purification", "","", "Random", ""),
        xlab="Scenario", ylab="Mean Delta Obs", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

# Var Delta
boxplot(simul$VarObs~simul$Conclusion*simul$Scenario, notch=F, outline=F,
        names=c("", "Positive", "", "", "Purification", "","", "Random", ""),
        xlab="Scenario", ylab="Var Delta Obs", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

boxplot(simul$Nmut~simul$Conclusion*simul$Scenario, notch=F, outline=F,
        names=c("", "Positive", "", "", "Purification", "","", "Random", ""),
        xlab="Scenario", ylab="Nb subsitution", col=col)
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

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
  obs[[TF]] <- data 
}

table(obs$Conclusion)/nrow(obs)*100
barplot(table(obs$Conclusion)/nrow(obs)*100, col=col, beside=F)

hist(obs$pval_NullPos, breaks=50)
hist(obs$pval_NullPos, breaks=100)
hist(obs$pval_PurifPos, breaks=100)

plot(log(obs$Nmut), col=col[obs$Conclusion], pch=16, cex=0.7)
boxplot(log(obs$Nmut)~obs$Conclusion)
