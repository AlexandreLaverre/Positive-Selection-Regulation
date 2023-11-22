simul <- read.csv("/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/simulation_results.csv", h=T)
simul$Conclusion <- as.factor(simul$Conclusion)
simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Null model"))

col=c("forestgreen", "deepskyblue3", "firebrick")
names(col) = c("Positive model", "Stabilizing model", "Null model")

pdf("/Users/alaverre/Documents/Detecting_positive_selection/results/MaxLikelihoodApproach/Comparisons_scenario_conclusion.pdf", width=7, height=8)
# Scenario ~ Conclusion
par(xpd = TRUE)
par(mfrow=c(2,1))
barplot(table(simul$Conclusion, simul$Scenario), col=col, las=1, xlab="Scenario", ylab="Nb simulation")
legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

text(x=2, y=1250, labels="Test conclusion")

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