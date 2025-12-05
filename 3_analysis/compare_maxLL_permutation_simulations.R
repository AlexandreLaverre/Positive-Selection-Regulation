library(corrplot)
library(RColorBrewer)
library(stringr)
library(dplyr)
par(xpd = TRUE)

sp="human"
TF="Wilson/HNF4A"
TFS=c("Wilson/HNF4A") #, "Wilson/FOXA1", "Wilson/HNF6", "Wilson/HNF4A", "Schmidt12/CTCF")
for (TF in TFS){
  print(TF)
  path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/", sp, "/", TF, "/")
  pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
  col=c("forestgreen", "deepskyblue3", "firebrick")
  names(col) = c("Positive model", "Stabilizing model", "Neutral model")
  
  new_col=c("forestgreen", "lightgreen", "deepskyblue3", "firebrick")
  names(new_col) =c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")
  
  col_permut=c("forestgreen", "lightgreen", "firebrick")
  names(col_permut) = c("Directional (+)", "Directional (-)", "Neutral")
  
  bin = "exact_ranked_50"
  type = "0.0Null_25Stab_25Pos_0.0Purif_addExtreme0"
  data = "beta" #"by_params_independent_SVM_without_backMut_noBin" #"by_deltas"by_params_independent_SVM_without_backMut_noBin
  epistasis = "NA" #ifelse(data=="by_params_independent_SVM_without_backMut_noBin", "without_epistasis", "with_epistasis")
  selection = c("stabilising", "positive", "neutral")
  maj = c("Stabilizing", "Positive", "Neutral")
  names(maj) <- selection
  maxSub = 150
  maxLen=1000
  
  ################################################################################
  # Read all results
  all_delta <- list()
  all_simul <- list()
  all_permut <- list()
  
  delta_total <- read.table(paste0(path, "/deltas/focal_all_possible_deltaSVM.txt"), h=T, sep="\t", quote="", fill=T, row.names = 1)
  #stats <- read.table(paste0(path, "/stats_simulation/simul_", data, ".txt"), h=T, sep="\t", row.names = 1)
  #stats$ratio <- stats$AlphaPos/stats$BetaPos
  #stats$classAlphaStab <- cut(stats$AlphaStab, breaks=quantile(stats$AlphaStab, probs=seq(0, 1, 0.1)), include.lowest = T)
  
  for (sel in selection){
    obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
    deltas <- read.table(paste0(path, "/deltas/simul_", data, "_", type, "_", sel, "_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
    row.names(deltas) <- deltas$seq_name
    all_delta[[paste0(sel, data)]] <- deltas
    
    # MaxLL
    simul <- read.csv(paste0(path, "/Tests/MLE_summary_simulated_", data, "_", type, "_", sel, "_", bin, "bins_threshold_0.01.csv"), h=T, row.names = 1)
    simul$Conclusion <- as.factor(simul$Conclusion)
    #simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Stabilizing model", "Neutral model"))
    #simul$Conclusion <- factor(simul$Conclusion, levels=c("Positive model", "Positive model Decreasing", "Positive model Increasing", "Stabilizing model", "Neutral model"))
    simul$Conclusion <- factor(simul$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))
  
    simul$ID <- rownames(simul)
    simul$Scenario <- maj[sel]
    simul$SVM <- deltas[rownames(simul), "SVM"]
    
    simul$deltaSVM <- deltas[row.names(simul),]$deltaSVM
    simul$Diffdelta <-simul$deltaSVM-simul$SumObs
    simul$length <- as.numeric(str_split_i(simul$ID, ":", 5))-as.numeric(str_split_i(simul$ID, ":", 4))
    simul$Prop.Sub <- simul$Nmut/simul$length
    
    all_simul[[paste0(sel, data)]] <- simul
    
    # Permutation Test
    permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_", data, "_", type, "_", sel, ".txt"), h=T, row.names = 1)
    permut <- permut[rownames(simul),]
    permut$signif <- as.factor(ifelse(permut$pval.high <= 0.01, "Pos Increase", ifelse(permut$pval.high >= 0.99, "Pos Decrease", "Neutral")))
    permut$signif.bin <- ifelse(permut$pval.high <= 0.01, 1, ifelse(permut$pval.high >= 0.99, -1, 0))
    permut$signif <- factor(permut$signif, levels=c("Pos Increase", "Pos Decrease", "Neutral"))
    
    #permut$signif <- as.factor(ifelse(permut$pval.high <= 0.01 | permut$pval.high >= 0.99, "Positive", "Neutral"))
    #permut$signif.bin <- ifelse(permut$pval.high <= 0.01 | permut$pval.high >= 0.99, 1, 0)
    #permut$signif <- factor(permut$signif, levels=c("Positive", "Neutral"))
    
    simul <- cbind(simul, permut[,c("SVM", "deltaSVM", "med.expected.deltaSVM", "mean.expected.deltaSVM", "pval.high", "signif", "signif.bin")])
    simul$LRT_direction <- as.factor(ifelse(simul$SumObs>simul$med.expected.deltaSVM, "Increasing", "Decreasing"))
    simul$Combine <- paste(simul$Conclusion, simul$LRT_direction)
    simul[which(simul$Conclusion=="Positive model"),]$Conclusion <- simul[which(simul$Conclusion=="Positive model"),]$Combine
  
    
    all_simul[[paste0(sel, data)]] <- simul
    
    #print(sel)
    #print(table(simul$Conclusion))
    a <- table(simul$Conclusion)
    
    if (sel == "positive"){
      message("TP: ", (a[1]+a[2])/10, "%")
    }else if(sel == "neutral"){
      message("FP: ", (a[1]+a[2])/10, "%")
    }
    #print(table(simul$signif))
  }
}


################################################################################
# Combine in one dataframe
rand <- all_simul[[paste0("neutral", data)]]
rownames(rand) <- NULL
pos <- all_simul[[paste0("positive", data)]]
rownames(pos) <- NULL
stab <- all_simul[[paste0("stabilising", data)]]
rownames(stab) <- NULL
simul <- rbind(rand, pos, stab)
simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Stabilizing", "Neutral"))

################################################################################
### Scenario to Conclusion
#pdf(paste0(pathFigure, "Simulation_comparisons_scenario_conclusion_", epistasis, ".pdf"), width=8, height=6)
par(xpd = TRUE)
par(mfrow=c(1,1))

# Total Scenario ~ Conclusion Max LL
barplot(table(simul$Conclusion, simul$Scenario), col=new_col, las=1, xlab="Scenario", 
        ylab="Nb simulation", cex.lab=1.4, cex.names=1.4, names=c("Directional", "Stabilising", "Random"))

legend("top", legend=c("Directional", "Stabilising", "Null"), fill=col, bty="n", ncol=3, inset = c(0, -0.17), cex=1.4)
text(x=2, y=1200, labels="Test Conclusion", cex=1.5)

# Total Scenario ~ Conclusion Test Permut
a <- table(simul$signif, simul$Scenario)
barplot(a, col=col_permut, las=1, xlab="Scenario", ylab="Nb simulation", cex.lab=1.5, cex.names=1.5)
legend("top", legend=names(col_permut), fill=col_permut, bty="n", ncol=2, inset = c(0, -0.25),  cex=1.5)
text(x=2, y=1300, labels="Permutations Test conclusion", cex=1.5)

# Overlap Tests
conclu.test.Max <- prop.table(table(simul$signif, simul$Conclusion))
barplot(conclu.test.Max, col=c("forestgreen","lightgreen", "firebrick"), las=1, ylab="Proportion", main="maxLL vs Permutation Test")

conclu.test <- prop.table(table(simul$Conclusion, simul$signif))
barplot(conclu.test, col=new_col, las=1, ylab="Proportion", main="Permutation vs maxLL Test")

#dev.off()

################################################################################
# Plot probabilities of substitutions
#pdf(paste0(pathFigure, "/Simulation_substitutions_probabilities_deltaSVM_no_epistasis.pdf"), width=8, height=7)
par(xpd = F)
par(mfrow=c(2,2))
for (sel in selection){
  proba_col = c("seq_id", paste("sub", 1:(maxLen*3), sep = ":"))
  proba <- read.table(paste0(path, "/substitution_probabilities/simul_", data, "_", sel, "_evolution.txt"), h=F, sep="\t", quote="", fill=T, col.names = proba_col, row.names = 1)
  
  seq = "chr2:33052075:33052325" # chr11:58762713:58762873
  proba_seq <- as.numeric(proba[seq,])
  proba_seq <- proba_seq[!is.na(proba_seq)]
  
  test_all <- as.numeric(delta_total[seq,])
  test_all <- test_all[!is.na(test_all)]
  obs_delta <- as.numeric(all_delta[[paste0(sel, data)]][seq,5:50])
  obs_delta <- obs_delta[which(!is.na(obs_delta))]
  
  plot(test_all, proba_seq, xlab="deltaSVM", ylab="Fixation probability", main=paste0(seq," - ", sel), cex=0.6)
  abline(v=obs_delta, col="red", lwd=0.5)
}

par(mfrow=c(2,1))
test_all <- t(delta_total[seq,])
test_all <- test_all[which(!is.na(test_all)),]
proba_seq <- proba[seq,]
proba_seq <- t(proba_seq[!is.na(proba_seq)])

delta_proba <- data.frame("deltaSVM"=test_all, "proba"=t(proba_seq))
delta_proba$direction <- sapply(strsplit(row.names(delta_proba), "\\."), function(x) x[2])
TT <- c(0.000469491339356585, 0.00219055617351051, 0.000541988047920129)
G <- c(0.00329670681456568, 0.000812337916957577, 0.000835488454479485)
C <- c(0.000815500132032239, 0.00328878175703288, 0.000836893883401556)
A <- c(0.000470545411048139, 0.000544369469149196, 0.00220437622457755)

delta_proba$proba.chr <- as.character(delta_proba$proba)
delta_proba$origin <- ifelse(delta_proba$proba.chr %in% TT, "T", ifelse(delta_proba$proba.chr %in% G, "G", ifelse(delta_proba$proba.chr %in% C, "C",  ifelse(delta_proba$proba.chr %in% A, "A", "Unknown"))))
delta_proba$sub <- paste0(delta_proba$origin, "-", delta_proba$direction)

delta_proba[which(delta_proba$origin %in% c("A", "T") & delta_proba$direction %in% c("A", "T")), "type"] <- "WW"
delta_proba[which(delta_proba$origin %in% c("A", "T") & delta_proba$direction %in% c("G", "C")), "type"] <- "WS"
delta_proba[which(delta_proba$origin %in% c("G", "C") & delta_proba$direction %in% c("G", "C")), "type"] <- "SS"
delta_proba[which(delta_proba$origin %in% c("G", "C") & delta_proba$direction %in% c("A", "T")), "type"] <- "SW"

boxplot(delta_proba$proba~delta_proba$sub, ylab="Fixation Probability", xlab="")
boxplot(delta_proba$proba~delta_proba$type, ylab="Fixation Probability", xlab="")

par(mfrow=c(2,2))
for (sel in selection){
  proba_col = c("seq_id", paste("sub", 1:(maxLen*3), sep = ":"))
  proba <- read.table(paste0(path, "/substitution_probabilities/simul_", data, "_", sel, "_evolution.txt"), h=F, sep="\t", quote="", fill=T, col.names = proba_col, row.names = 1)
  Nb_proba <- apply(proba, 1, function(x) length(x[!is.na(x) & as.numeric(x)>0])) 
  Nb_possible <- apply(proba, 1, function(x) length(x[!is.na(x)]))
  hist(Nb_proba, breaks=50, xlab="Nb Sub with p(fix)>0", main=sel)
  mtext(paste0("min NbSub=", min(Nb_proba)))
  hist(Nb_proba/Nb_possible, breaks=50, xlab="Seq proportion with p(fix)>0", main=sel)
}

dev.off()

################################################################################

# simul$Conclusion <- NULL  by_params_epistasis_without_backMut
# CorMat <- cor(simul)
# corrplot(CorMat, type="upper", order="hclust", main=type, col=brewer.pal(n=8, name="RdYlBu"))

################################################################################
# Effect of features
for (sel in selection){
  scenar <- all_simul[[paste0(sel, data)]]
  scenar$classSVM <- cut(scenar$SVM, breaks=quantile(scenar$SVM, probs = seq(0, 1, 0.1)), include.lowest = T)
  scenar$classdeltaSVM <- cut(scenar$deltaSVM, breaks=quantile(scenar$deltaSVM, probs = seq(0, 1, 0.1)), include.lowest = T)
  scenar$Model.bin <- ifelse(scenar$Conclusion == paste(maj[[sel]], "model"), 1, 0)
  scenar$positiveDelta <- ifelse(scenar$deltaSVM > 0, 1, 0)
  
  scenar$positive <- ifelse(simul$Conclusion %in% c("Directional (+)", "Directional (-)"), 1, 0)
  scenar$neutral <- ifelse(simul$Conclusion == "Neutral model", 1, 0)
  simul$stabilising <- ifelse(simul$Conclusion == "Stabilizing", 1, 0)
  # Effect of substitutions number
  par(mfrow=c(1,2))
  a <- table(scenar$Nmut, scenar$signif.bin)
  plot(a[,2]/(a[,2]+a[,1]), las=1, type="b", xlab="Nb substitutions", ylab="Proportion signif", main=paste("Permut Test", sel))
  
  a <- table(scenar$Nmut, scenar$Model.bin)
  barplot(a[,2]/(a[,2]+a[,1]), las=1,  type="b", xlab="Nb substitutions", ylab="Proportion signif", main=paste("LL Ratio Test", sel))
  
  # Effect of SVM quantile
  par(mfrow=c(1,2))
  a <- table(scenar$classSVM,scenar$signif.bin)
  barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="SVM quantile", ylab="Proportion signif", main="Permut Test")
  
  a <- table(scenar$classSVM, scenar$Model.bin)
  barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="SVM quantile", ylab="Proportion signif", main="LL Ratio Test")
  
  # Effect of deltaSVM quantile
  par(mfrow=c(1,2))
  a <- table(scenar$classdeltaSVM,scenar$signif.bin)
  barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="deltaSVM quantile", ylab="Proportion signif", main="Permut Test")
  
  a <- table(scenar$classdeltaSVM, scenar$Model.bin)
  barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="deltaSVM quantile", ylab="Proportion signif", main="LL Ratio Test")
  
  # Plot correlation Apha Simul versus Alpha Estimated
  stats <- stats[rownames(scenar),]
  scenar$ratio <- scenar$AlphaPos/scenar$BetaPos
  
  if (sel == "positive"){
    # Correlation between ratio simulated and estimated
    seq_alpha_pos <- rownames(stats[which(stats$AlphaPos > 1),])
    seq_alpha_neg <- rownames(stats[which(stats$AlphaPos == 1),])
    boxplot(log(scenar[seq_alpha_pos,]$ratio)~log(stats[seq_alpha_pos,]$ratio), outline=F, xlab="Alpha Simul", ylab="Alpha Estim")
    boxplot(log(scenar[seq_alpha_neg,]$ratio)~log(stats[seq_alpha_neg,]$ratio), outline=F, xlab="Alpha Simul", ylab="Alpha Estim")
    plot(log(stats$ratio), log(scenar$ratio), xlab="log(Alpha/Beta) Simulated", ylab="log(Alpha/Beta) Estimated")
    
    # Effect of Alpha/Beta value on positive proportion
    a <- table(stats$ratio, scenar$Model.bin)
    barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="AlphaPos/BetaPos", ylab="Proportion signif", main="LL Ratio Test")
  }
  
  if (sel == "stabilising"){
    # Correlation of simulated and estimated Alpha
    plot(log(stats$AlphaStab), log(scenar$AlphaPurif))
    boxplot(log(scenar$AlphaPurif)~stats$classAlphaStab, outline=F)
    
    # Effect of Alpha on stabilising proportion
    a <- table(stats$classAlphaStab, scenar$Model.bin)
    barplot(a[,2]/(a[,2]+a[,1]), las=1, xlab="class AlphaStab", ylab="Proportion signif", main="LL Ratio Test")
  }
}


###### Check params effect
scenar$AlphaPosSimul <- stats[row.names(scenar),]$ratio

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Impact SVM
par(mfrow=c(1,2))
a <- table(scenar$classSVM, scenar$signif.bin)/nrow(scenar)*10
barplot(t(a), col=c("firebrick", "forestgreen"),
        las=1, xlab="SVM quantile", ylab="Proportion", main="Permut Test")

a <- table(scenar$classdeltaSVM, scenar$signif.bin)/nrow(scenar)*10
barplot(t(a), col=c("firebrick", "forestgreen"),
        las=1, xlab="deltaSVM quantile", ylab="Proportion", main="Permut Test")


a <- table(scenar$classSVM, scenar$Conclusion)/nrow(scenar)*10
barplot(t(a), col=c("forestgreen", "deepskyblue3", "firebrick"), las=1, xlab="SVM quantile", 
        ylab="Proportion", main="LL Ratio Test")

a <- table(scenar$classdeltaSVM, scenar$Conclusion)/nrow(scenar)*10
barplot(t(a), col=c("forestgreen", "deepskyblue3", "firebrick"),
        las=1, xlab="deltaSVM quantile", ylab="Proportion", main="LL Ratio Test")

################################################################################ 

pdf(paste0(path, "figures/Epistatic_effect_simulations_", data, ".pdf"))
par(mfrow=c(2,2))
hist(simul$Diffdelta/(abs(simul$SumObs)+abs(simul$deltaSVM))/2, breaks=100, main="All simulations", xlab="Difference Delta / mean(SumDelta, AllDelta)")
hist(pos$Diffdelta/(abs(pos$SumObs)+abs(pos$deltaSVM))/2, breaks=100, main="Positive simulations", xlab="Difference Delta / mean(SumDelta, AllDelta)")
hist(rand$Diffdelta/(abs(rand$SumObs)+abs(rand$deltaSVM))/2, breaks=100, main="Random simulations", xlab="Difference Delta / mean(SumDelta, AllDelta)")
hist(stab$Diffdelta/(abs(stab$SumObs)+abs(stab$deltaSVM))/2, breaks=100, main="Stabilising simulations", xlab="Difference Delta / mean(SumDelta, AllDelta)")
dev.off()

boxplot(simul$Diffdelta~simul$Scenario, notch=T, outline=F, main=data,
        xlab="Scenario", ylab="Diff delta", col=col)

################################################################################
par(mfrow=c(2,2))
# Mean Delta
boxplot(simul$MeanObs~simul$Conclusion*simul$Scenario, notch=T, outline=F,
        names=c("", "Positive", "", "", "Stabilizing", "","", "Random", ""),
        xlab="Scenario", ylab="Mean Delta Obs", col=col)
legend("top", legend=names(col), fill=col, bty="n", cex=0.8, ncol=3, inset = c(0, -0.15))

# Var Delta
boxplot(simul$VarObs~simul$Conclusion*simul$Scenario, notch=F, outline=F,
        names=c("", "Positive", "", "", "Stabilizing", "","", "Random", ""),
        xlab="Scenario", ylab="Var Delta Obs", col=col)
#legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))

boxplot(simul$Nmut~simul$Conclusion*simul$Scenario, notch=F, outline=F,
        names=c("", "Positive", "", "", "Stabilizing", "","", "Random", ""),
        xlab="Scenario", ylab="Nb subsitution", col=col)

boxplot(simul$Nmut~simul$signif*simul$Scenario, notch=T, outline=F,
        xlab="Scenario", ylab="Nb subsitution", col=col_permut)

#legend("top", legend=names(col), fill=col, bty="n", ncol=3, inset = c(0, -0.25))
abline(v=c(3.5,6.5), lwd=0.8)
dev.off()

################################################################################
## Epistatic effect
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



################################################################################
# Delta per position
delta_per_pos <- as.data.frame(t(delta_total))
delta_per_pos$pos <- str_split_i(row.names(delta_per_pos), "\\.", 1)
result <- tapply(delta_per_pos, delta_per_pos$pos, median)

med <- list()
conf_low <- list()
conf_sup <- list()
for (pos in unique(delta_per_pos$pos)){
  a <- boxplot(unlist(delta_per_pos[which(delta_per_pos$pos == pos),-1001]), na.rm=T, plot=F)
  med[[pos]] <- a$stats[3,]
  conf_low[[pos]] <- a$conf[1,]
  conf_sup[[pos]] <- a$conf[2,]
}

pdf(paste0(pathFigure, "/deltaSVM_along_positions.pdf"))
par(xpd = F)
df <- data.frame("median"=unlist(med), "conf_low"=unlist(conf_low), "conf_high"=unlist(conf_sup))
df <- df[0:847,]
plot(df$median, cex=0.5, type="l", xlab="Position on sequence", ylab="deltaSVM (median)", main="All possible deltaSVM on 1000 human CEBPA TFBS")

polygon(c(seq_along(df$median), rev(seq_along(df$median))), c(df$conf_low, rev(df$conf_high)), col = rgb(1, 0, 0, 0.4), border=NA)
dev.off()

################################################################################
alpha= 6.163961e-06
beta=1.8877228
par(mfrow=c(1,1))

### 
p = seq(0, 1, length=1000)
#create plot of Beta distribution with shape parameters 2 and 10
plot(p, dbeta(p, alpha, beta), type='l', ylab="Relative Fixation Probability", xlab=expression(Delta*" SVM"), 
     las=1, col="forestgreen", lwd=3, cex.axis=1, cex.lab=1.2, main="Beta distributions", axes=T)


lines(p, dbeta(p, 40, 2), col="navy", lwd=3)
lines(p, rep(5, length(p)), col="firebrick", lwd=3)
axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(-5, -2.5, 0, 2.5, 5))
axis(2, at=c(0, 2.5, 5, 7.5, 10), labels=c(0, 0.5, 1, 1.5, 2), las=2)
legend("topleft", legend=c("Beta=2", "Beta=1"), lty=c(1,1,1), lwd=3,
       col=c("navy", "forestgreen"), cex=1.2)

