# Figure simulation
library(RColorBrewer)
par(xpd = TRUE)

sp="drosophila"
TF="Ni12/CTCF"
path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/", sp, "/", TF, "/")
pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

new_col=c("#238b45", "#bae4b3", "#4393C3", "#D6604D")
names(new_col) =c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")

col_permut=c("#238b45", "#bae4b3", "#D6604D")
names(col_permut) =c("Directional (+)", "Directional (-)", "Neutral model")

bin = "exact_ranked_50"
selection = c("stabilising", "positive", "neutral")
maj = c("Stabilizing", "Positive", "Random")
names(maj) <- selection
maxSub = 150
maxLen = 1000

params_pos=c(2, 5, 10, 25, 50, 100) # pos
params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
params_purif=c("0.0", 0.1, 0.2, 0.5, "1.0", "2.0", "3.0", "5.0", "10.0")
addExtreme = "0.0Purif_addExtreme0_"

#params_pos=c(25) # pos
#params_null=c("0.0")
#params_purif=c("0.0")
#addExtreme = "_0.0Purif_addExtreme1_"

get_all_simul <- function(params, type="pos", sel="positive"){
  all_delta <- list()
  all_simul <- list()
  
  for (param in params){
    if (type == "null"){
      data = paste0("beta_", param, "Null_25Stab_25Pos_")
    }else if (type == "pos"){
      data = paste0("beta_0.0Null_25Stab_", param, "Pos_")
    }else if (type == "purif"){
      data = paste0("beta_0.0Null_25Stab_25Pos_", param, "Purif_")
    }
    
    obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
    deltas <- read.table(paste0(path, "/deltas/simul_", data, addExtreme, sel, "_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
    row.names(deltas) <- deltas$seq_name
    all_delta[[paste0(sel, data)]] <- deltas
    
    simul <- read.csv(paste0(path, "/Tests/MLE_summary_simulated_", data, addExtreme, sel, "_", bin, "bins_threshold_0.01.csv"), h=T, row.names = 1)
    simul$Conclusion <- as.factor(simul$Conclusion)
    simul$Conclusion <- factor(simul$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))
    simul$ID <- rownames(simul)
    simul$Scenario <- maj[sel]
    simul$param <- param
    
    simul$SVM <- deltas[rownames(simul), "SVM"]
    simul$deltaSVM <- deltas[row.names(simul),]$deltaSVM
    
    # Permut
    permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_", data, addExtreme, sel, ".txt"), h=T, row.names = 1)
    permut <- permut[rownames(simul),]
    permut$pval.two.tailed <- apply(permut, 1, function(x) 2*min(as.numeric(x["pval.high"]), 1-as.numeric(x["pval.high"])))
    permut$signif <- as.factor(ifelse(permut$pval.two.tailed > 0.01, "Neutral", ifelse(permut$deltaSVM >0, "Directional (+)", "Directional (-)")))
    permut$signif.bin <- ifelse(permut$pval.two.tailed <= 0.01, 1, 0)
    permut$signif <- factor(permut$signif, levels=c("Directional (+)", "Directional (-)", "Neutral"))
    
    simul <- cbind(simul, permut[,c("pval.two.tailed", "signif", "signif.bin")])
    all_simul[[paste0(param)]] <- simul
  }
  
  all <- do.call(rbind, all_simul)
  
  return(all)
}


#pdf(paste0(pathFigure, "Positive_Rate_simulations_", sel, ".pdf"))
par(mfrow=c(2,2))

sel="positive"
all = get_all_simul(params_null, type="null", sel=sel)

# VarObs
max = ifelse(sel == "neutral", 1, 0.1)
ylab = ifelse(sel == "positive",   "Proportion True Positive", "Proportion False Positive")
all$classVarObs <- cut(all$deltaSVM, breaks=quantile(all$VarObs, probs = seq(0, 1, 0.05)), include.lowest = T)
prop_TP_VarObs = prop.table(table(all$Conclusion, all$classVarObs), margin=2)
plot(prop_TP_VarObs[1,]+prop_TP_VarObs[2,], xlab="quantile VarObs", ylab= ylab, type="b", xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_VarObs[2,]), labels=names(prop_TP_VarObs[2,]))

prop_TP_VarObs = prop.table(table(all$signif, all$classVarObs), margin=2)
lines(prop_TP_VarObs[2,]+prop_TP_VarObs[1,], col="red", type="b")
legend("top", legend=c("MLL", "Permutations"), col=c("black", "red"), bty="n", lty=1)


# deltaSVM
max = ifelse(sel == "neutral", 1, 0.1)
ylab = ifelse(sel == "positive",   "Proportion True Positive", "Proportion False Positive")
all$classdeltaSVM <- cut(all$deltaSVM, breaks=quantile(all$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)
prop_TP_delta = prop.table(table(all$Conclusion, all$classdeltaSVM), margin=2)
plot(prop_TP_delta[1,]+prop_TP_delta[2,], xlab="quantile deltaSVM", ylab= ylab, type="b", xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))

prop_TP_delta = prop.table(table(all$signif, all$classdeltaSVM), margin=2)
lines(prop_TP_delta[2,]+prop_TP_delta[1,], col="red", type="b")
legend("top", legend=c("MLL", "Permutations"), col=c("black", "red"), bty="n", lty=1)

# SVM
max = ifelse(sel == "neutral", 0.05, 0.1)
ylab = ifelse(sel == "neutral", "Proportion False Positive",  "Proportion True Positive")
all$classSVM <- cut(all$SVM, breaks=quantile(all$SVM, probs = seq(0, 1, 0.05)), include.lowest = T)
prop_TP_SVM = prop.table(table(all$Conclusion, all$classSVM), margin=2)
plot(prop_TP_SVM[1,]+prop_TP_SVM[2,], xlab="quantile SVM", ylab= ylab, type="b", xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_SVM[2,]), labels=names(prop_TP_SVM[2,]))

prop_TP_SVM = prop.table(table(all$signif, all$classSVM), margin=2)
lines(prop_TP_SVM[2,]+prop_TP_SVM[1,], col="red", type="b")
legend("top", legend=c("MLL", "Permutations"), col=c("black", "red"), bty="n", lty=1)


# Substitutions 
max = ifelse(sel == "neutral", 0.1, 0.1)
prop_TP_sub = prop.table(table(all$Conclusion, all$Nmut), margin=2)
plot(prop_TP_sub[2,]+prop_TP_sub[1,], xlab="Nb Substitution", ylab= ylab, type="b", xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_sub[2,]), labels=names(prop_TP_sub[2,]))
prop_TP_sub = prop.table(table(all$signif, all$Nmut), margin=2)
lines(prop_TP_sub[2,]+prop_TP_sub[1,], col="red", type="b")

# Positive Strenght
#all$param <- as.factor(all$param)
max = ifelse(sel == "neutral", 0.02, 1)
prop_TP_delta = prop.table(table(all$Conclusion, all$param), margin=2)
plot(prop_TP_delta[1,]+prop_TP_delta[2,], xlab="Proportion of Null", ylab= ylab, type="b",  xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))

prop_TP_delta = prop.table(table(all$signif, all$param), margin=2)
lines(prop_TP_delta[2,]+prop_TP_delta[1,], col="red", type="b")

# Directional Selection Strength
all = get_all_simul(params_pos, type="pos", sel=sel)
max = ifelse(sel == "neutral", 0.02, 1)
prop_TP_delta = prop.table(table(all$Conclusion, all$param), margin=2)
plot(prop_TP_delta[1,]+prop_TP_delta[2,], xlab="Directional Selection Strength", ylab= ylab, type="b",  xaxt="n", ylim=c(0,max))
axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))

prop_TP_delta_signif = prop.table(table(all$signif, all$param), margin=2)
lines(prop_TP_delta[2,]+prop_TP_delta[1,], col="red", type="b")

dev.off()

all = get_all_simul(params_purif, type="purif", sel="neutral")
all$param <- factor(all$param, levels=params_purif)
prop_TP_delta = prop.table(table(all$Conclusion, all$param), margin=2)
plot(prop_TP_delta[1,]+prop_TP_delta[2,], xlab="Scale Purif", ylab= ylab, type="b",  xaxt="n", ylim=c(0,0.55))
axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))

prop_TP_delta = prop.table(table(all$signif, all$param), margin=2)
lines(prop_TP_delta[2,]+prop_TP_delta[1,], col="red", type="b")

# Plot Proportion False Positive per param and per substitutions
pdf(paste0(pathFigure, "FP_with_purifying_selection.pdf"))
par(mfrow=c(2,2))
for (i in params_purif){
  subdata = all[which(all$param == i),]
  cor <- cor.test(subdata$deltaSVM,subdata$Nmut)$estimate
  boxplot(subdata$deltaSVM~subdata$Nmut, main="Neutral", xlab="Nb Substitution", ylab= "deltaSVM")
  mtext(paste("R=", round(cor, 2)), side=3)

  prop_TP_sub = prop.table(table(subdata$Conclusion, subdata$Nmut), margin=2)
  plot(prop_TP_sub[2,]+prop_TP_sub[1,], xlab="Nb Substitution", ylab= "Proportion False Positive", type="b", xaxt="n")
  axis(1,at=1:length(prop_TP_sub[2,]), labels=names(prop_TP_sub[2,]))
  prop_TP_sub = prop.table(table(subdata$signif, subdata$Nmut), margin=2)
  lines(prop_TP_sub[2,]+prop_TP_sub[1,], col="red", type="b")
  legend("topleft", legend=c("MLL", "Permut"), col=c("black", "red"), bty="n", lty=1)

}
dev.off()


