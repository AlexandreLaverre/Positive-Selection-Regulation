# Figure simulation

library(corrplot)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
par(xpd = TRUE)

sp="human"
TF="Wilson/HNF4A"
peakType="NarrowPeaks"
path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/", peakType, "/", sp, "/", TF, "/")
pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

new_col=c("#238b45", "#bae4b3", "#4393C3", "#D6604D")
names(new_col) =c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")

col_permut=c("#238b45", "#bae4b3", "#D6604D")
names(col_permut) =c("Directional (+)", "Directional (-)", "Neutral model")

bin = "exact_absolute_50"
data = "beta_0.1Null_25Stab_25Pos"
selection = c("stabilising", "positive", "neutral")
maj = c("Stabilizing", "Positive", "Random")
names(maj) <- selection
maxSub = 150
maxLen = 1000

################################################################################
# Read all results
all_delta <- list()
all_simul <- list()
delta_total <- read.table(paste0(path, "/deltas/focal_all_possible_deltaSVM.txt"), h=T, sep="\t", quote="", fill=T, row.names = 1)

for (sel in selection){
  obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
  
  #if (grepl("positive", sel)){sel_d=str_split_i(sel, "_", 2)}else{sel_d=sel}
  sel_d=sel
  deltas <- read.table(paste0(path, "/deltas/simul_", data, "_", sel_d, "_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
  row.names(deltas) <- deltas$seq_name
  all_delta[[paste0(sel, data)]] <- deltas
  
  # MaxLL
  simul <- read.csv(paste0(path, "/Tests/MLE_summary_simulated_", data, "_", sel, "_", bin, "bins_threshold_0.01.csv"), h=T, row.names = 1)
  simul$New_Conclusion <- simul$Conclusion 
  simul$New_Conclusion <- ifelse(simul$p_value_null_pos<0.1 & simul$p_value_null_purif>0.1 & simul$p_value_purif_pos>0.1, "Low_Dir", simul$Conclusion)
  simul$New_Conclusion <- as.factor(simul$New_Conclusion)
  simul$Conclusion <- factor(simul$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model", "Low_Dir"))
  
  
  simul$Conclusion <- as.factor(simul$Conclusion)
  simul$Conclusion <- factor(simul$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model", "Low_Dir"))


  simul$ID <- rownames(simul)
  simul$Scenario <- maj[sel]
  simul$SVM <- deltas[rownames(simul), "SVM"]
  
  simul$deltaSVM <- deltas[row.names(simul),]$deltaSVM
  simul$Diffdelta <-simul$deltaSVM-simul$SumObs
  simul$length <- as.numeric(str_split_i(simul$ID, ":", 5))-as.numeric(str_split_i(simul$ID, ":", 4))
  simul$Prop.Sub <- simul$Nmut/simul$length
  
  # Sequences
  sequences <- readDNAStringSet(paste0(path, "/sequences/simulated_sequences_by_", data, "_", sel, "_evolution.fa"), format = "fasta")
  
  # Permutation Test
  
  permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_", data, "_", sel, ".txt"), h=T, row.names = 1)
  #permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_by_params_independent_SVM_without_backMut_noBin_", sel, ".txt"), h=T, row.names = 1)
  
  permut <- permut[rownames(simul),]
  permut$signif <- as.factor(ifelse(permut$pval.high <= 0.01, "Directional (+)", ifelse(permut$pval.high >= 0.99, "Directional (-)", "Neutral")))
  permut$signif.bin <- ifelse(permut$pval.high <= 0.01, 1, ifelse(permut$pval.high >= 0.99, 1, 0))
  permut$signif <- factor(permut$signif, levels=c("Directional (+)", "Directional (-)", "Neutral"))
  
  simul <- cbind(simul, permut[,c("pval.high", "signif", "signif.bin")])
  all_simul[[paste0(sel, data)]] <- simul
  
  a <- table(simul$Conclusion)

  if (grepl("positive", sel)){
    message(sel)
    message("TP: ", (a[1]+a[2])/50, "%")
  }else if(sel == "neutral"){
    message("FP: ", (a[1]+a[2])/50, "%")
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
simul$Scenario <- factor(simul$Scenario, levels=c("Positive", "Positive -", "Stabilizing", "Random"))
simul[which(simul$SumObs<0 & simul$Scenario=="Positive"),]$Scenario <- "Positive -"

table(simul$Conclusion, simul$Scenario)
table(simul$New_Conclusion, simul$Scenario)

################################################################################
### Scenario to Conclusion
simp_TF=strsplit(TF, "/")[[1]][2]
pdf(paste0(pathFigure, "Simulation_", sp, "_", simp_TF, "_",  data, "_", bin, "", ".pdf"), width=8, height=6)
par(xpd = TRUE)
par(mfrow=c(2,2), mar = c(4, 4, 4, 1))

# Max LH
barplot(prop.table(table(simul$Conclusion, simul$Scenario), margin=2), col=new_col, las=1, xlab="Scenario",
        ylab="Proportion", cex.lab=1, cex.names=1, names=c("Directional (+)", "Directional (-)", "Stabilising", "Random"))

legend("top", legend=c("Directional (+)", "Directional (-)", "Stabilising", "Random"), 
       fill=new_col, bty="n", ncol=2, inset = c(0, -0.35), cex=1)
text(x=2, y=1200, labels="Test Conclusion", cex=1.5)

# Test Permut
barplot(prop.table(table(simul$signif, simul$Scenario), margin=2), col=col_permut, names=c("Directional (+)", "Directional (-)", "Stabilising", "Random"),
        las=1, xlab="Scenario", ylab="Nb simulation", main=paste(sp, TF, "Permutations Test"))

# Overlap Tests
conclu.test.Max <- prop.table(table(simul$signif, simul$Conclusion))
barplot(conclu.test.Max, col=c("#238b45", "#bae4b3", "#D6604D"), las=1, ylab="Proportion", main="maxLL vs Permutation Test")

conclu.test <- prop.table(table(simul$Conclusion, simul$signif))
barplot(conclu.test, col=new_col, las=1, ylab="Proportion", main="Permutation vs maxLL Test")

################################################################################
for (sel in selection){
  par(mfrow=c(1,2))
  simul <- all_simul[[paste0(sel, data)]]
  simul$positive <- ifelse(simul$Conclusion %in% c("Directional (+)", "Directional (-)"), 1, 0)
  simul$neutral <- ifelse(simul$Conclusion == "Neutral model", 1, 0)
  simul$stabilising <- ifelse(simul$Conclusion == "Stabilizing", 1, 0)

  simul$permut_positive <- simul$signif.bin
  simul$permut_neutral <-  ifelse(simul$signif == "Neutral", 1, 0)

  ##### Proportion TP vs nb substitutions 
  #if (grepl("positive", sel)){sel_d=str_split_i(sel, "_", 2)}else{sel_d=sel}
  sel_d=sel
  sel_perm=paste0("permut_", sel)
  #simul$classSub <- cut(simul$Nmut, breaks=quantile(simul$Nmut, probs = seq(0, 1, 0.1)), include.lowest = T)
  # Max LL
  prop_TP_sub = prop.table(table(simul[,sel_d], simul$Nmut), margin=2)
  plot(prop_TP_sub[2,], xlab="Nb Substitution", ylab= "Proportion True Positive", type="b", main=sel, xaxt="n")
  axis(1,at=1:length(prop_TP_sub[2,]), labels=names(prop_TP_sub[2,]))
  
  # Permut
  if (sel!="stabilising"){
   prop_TP_sub = prop.table(table(simul[,sel_perm], simul$Nmut), margin=2)
   lines(prop_TP_sub[2,], col="red", type="b")
    #plot(prop_TP_sub[2,], xlab="Nb Substitution", ylab= "Proportion True Positive", type="b", main=sel_perm, xaxt="n")
    #axis(1,at=1:length(prop_TP_sub[2,]), labels=names(prop_TP_sub[2,]))
  }
  
  # Proportion TP vs deltaSVM
  simul$classdeltaSVM <- cut(simul$deltaSVM, breaks=quantile(simul$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)
  
  prop_TP_delta = prop.table(table(simul[,sel_d], simul$classdeltaSVM), margin=2)
  plot(prop_TP_delta[2,], xlab="quantile deltaSVM", ylab= "Proportion True Positive", type="b", main=sel, xaxt="n", ylim=c(0,1))
  axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))

  # Permut
  if (sel!="stabilising"){
    simul$classdeltaSVM <- cut(simul$deltaSVM, breaks=quantile(simul$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)
    
    prop_TP_delta = prop.table(table(simul[,sel_perm], simul$classdeltaSVM), margin=2)
    lines(prop_TP_delta[2,], col="red", type="b")
    #plot(prop_TP_delta[2,], xlab="quantile deltaSVM", ylab= "Proportion True Positive", type="b", main=sel_perm, xaxt="n", ylim=c(0,1))
    #axis(1,at=1:length(prop_TP_delta[2,]), labels=names(prop_TP_delta[2,]))
  }

  ###### Proportion TP vs SVM
  #simul$classSVM <- cut(simul$SVM, breaks=quantile(simul$SVM, probs = seq(0, 1, 0.05)), include.lowest = T)
  #prop_TP_SVM = prop.table(table(simul[,sel_d], simul$classSVM), margin=2)
  #plot(rev(prop_TP_SVM[2,]), xlab="quantile SVM", ylab= "Proportion True Positive", type="b", main=sel, xaxt="n", ylim=c(0.75,1))
  #axis(1,at=1:length(prop_TP_SVM[2,]), labels=rev(names(prop_TP_SVM[2,])))
  
}

### Heatmaps
# Substitution to deltaSVM
simul <- rbind(rand, pos, stab)
simul$SumObs_class <- cut(simul$SumObs,  breaks=quantile(simul$SumObs, probs = seq(0, 1, 0.05)), include.lowest = T)
simul$Sub_class <- cut(simul$Prop.Sub,  breaks=quantile(simul$Prop.Sub, probs = seq(0, 1, 0.05)), include.lowest = T)

col <- colorRampPalette(brewer.pal(8, "Spectral"))(25)

simul$Pos.bin <- ifelse(grepl("Directional", simul$Conclusion), 1, 0)
simul$Stab.bin <- ifelse(simul$Conclusion=="Stabilizing", 1, 0)
simul$Rand.bin <- ifelse(simul$Conclusion=="Neutral model", 1, 0)


# Heatmap 
# Calculate proportions
proportion <- as.data.frame(scale(t(table(simul$SumObs_class, simul$Sub_class))))
colnames(proportion) <- c("SumObs_class", "Sub_class", "Frequency")
proportion_stab <- simul %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Stab.bin))
proportion_pos <- simul %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Pos.bin))
proportion_rand <- simul %>% group_by(SumObs_class, Sub_class) %>% summarise(Proportion = mean(Rand.bin))

#pdf(paste0(path, "figures/All_peaks_proportion_delta_sub_heatmap.pdf"), width = 10)
par(mfrow=c(1,1))
#ggplot(proportion, aes(x = Sub_class, y = SumObs_class, fill = Frequency)) +
#  geom_tile() +   scale_fill_gradientn(colors=rev(col)) + theme(legend.position = "none") + 
#  labs(title = "Scaled frequency of peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_pos, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "forestgreen") + 
  labs(title = "Proportion of Positive peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_stab, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "deepskyblue3") + 
  labs(title = "Proportion of Stabilising peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

ggplot(proportion_rand, aes(x = SumObs_class, y = Sub_class, fill = Proportion)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "firebrick") + 
  labs(title = "Proportion of Neutral peaks", x = "Sum deltaSVM per substitution", y = "Substitution per bp")

#dev.off()



par(mfrow=c(2,2))
hist(rand$p_value_null_purif, breaks=100, main="Null vs Stabilising", xlab="p-values")
hist(rand$p_value_null_pos, breaks=100, main="Null vs Directional", xlab="p-values")
hist(rand$p_value_purif_pos, breaks=100, main="Stabilising vs Directional", xlab="p-values")

dev.off()
