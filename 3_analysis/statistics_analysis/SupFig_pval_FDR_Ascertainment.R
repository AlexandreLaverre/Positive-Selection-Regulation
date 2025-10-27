library(ggplot2)
library(gridExtra)
library(grid)
library(qvalue)
library(reshape2)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
pathFigures <- paste0(path, "../../final_figures/")
sp = "human"
TF = "Wilson/CEBPA" #c("Wilson/CEBPA", "Wilson/HNF6", "Wilson/HNF4A", "Wilson/FOXA1") #, "Schmidt12/CTCF")
method = ifelse(sp=="human", "exact_ranked_ancestral", "exact")
cutoff = 2
treshold = ifelse(sp=="drosophila", 0.05, 0.1)
subMax = ifelse(sp=="drosophila", 50, 13)
leg_x = ifelse(sp=="drosophila", 0.4, 0.6)
bias <- list()

sweep_cutoffs <- function(maxLL, pvec_name="p_value_purif_pos", Ncol="N_substitution", cutoffs=c(1:10), thr=treshold) {
  pvec <- maxLL[[pvec_name]]
  pvec_permut <- maxLL$pvalue_perm
  Nsub <- maxLL[[Ncol]]
  out <- data.frame(cutoff = cutoffs, n_tests = NA, n_sig = NA, perm=NA, perm_FDR=NA, n_sigFDR=NA, pi0 = NA, perc = NA, perc_perm = NA, perc_perm_FDR = NA, percFDR = NA)
  for (i in seq_along(cutoffs)) {
    k <- cutoffs[i]
    keep <- which(Nsub >= k)
    out$n_tests[i] <- length(keep)
    if (length(keep) == 0) { out$n_sig[i] <- 0; out$pi0[i] <- NA; next }
    qval <- qvalue(pvec[keep])
    fdr <- p.adjust(pvec[keep], method = "BH")
    fdr_permut <- p.adjust(pvec_permut[keep], method = "BH")
    out$perm[i] <- sum(pvec_permut[keep] <= 0.01, na.rm = TRUE)
    out$perm_FDR[i] <- sum(fdr_permut <= thr, na.rm = TRUE)
    out$n_sig[i] <- sum(qval$q <= thr, na.rm = TRUE)
    out$n_sigFDR[i] <- sum(fdr <= thr, na.rm = TRUE)
    out$pi0[i] <- round(qval$pi0, digits=2)
    out$perc_perm[i] <- round(100*out$perm[i] / out$n_tests[i], digits=2)
    out$perc_perm_FDR[i] <- round(100*out$perm_FDR[i] / out$n_tests[i], digits=2)
    out$perc[i] <- round(100*out$n_sig[i] / out$n_tests[i], digits=2)
    out$percFDR[i] <- round(100*out$n_sigFDR[i] / out$n_tests[i], digits=2)
  }
  return(out)
}


print(TF)
file = paste0(path, "positive_selection/NarrowPeaks/", sp, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
fileMaxLL = paste0(path, "positive_selection/NarrowPeaks/", sp, "/", TF, "/Tests/MLE_summary_", method, ".csv")

################################################################################
maxLL <- read.csv(fileMaxLL, h=T, row.names = 1)

rand <- read.table(file, h=T)
rownames(rand) <- rand$ID
rand <- rand[match(rownames(maxLL), rand$ID),]
maxLL$pvalue_perm <- rand$pval.high

sweep <- sweep_cutoffs(maxLL, pvec_name="p_value_purif_pos", Ncol="Nmut", cutoffs=2:10, thr=treshold)
# Reshape for plotting
df_plot <- reshape2::melt(sweep, id.vars = "cutoff", measure.vars = c("perm", "perm_FDR", "n_sigFDR"))
df_plot$variable <- factor(df_plot$variable, levels = c("perm", "perm_FDR", "n_sigFDR"))

# Create diagnostic plot
p1 <- ggplot(data.frame(p = maxLL$p_value_null_purif), aes(x = p)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Null vs Stabilising", x = "p-values", y = "Count") +
  theme_minimal(base_size = 11)

p2 <- ggplot(data.frame(p = maxLL$p_value_purif_pos), aes(x = p)) +
  geom_histogram(bins = 50, fill = "salmon", color = "black") +
  labs(title = "Stabilising vs Directional", x = "p-values", y = "Count") +
  theme_minimal(base_size = 11)

p3 <- ggplot(data.frame(p = maxLL$p_value_null_pos), aes(x = p)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "black") +
  labs(title = "Null vs Directional", x = "p-values", y = "Count") +
  theme_minimal(base_size = 11)

p_sweep <- ggplot(df_plot, aes(x = cutoff, y = value, color = variable, group = variable)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("perm" = "black", "perm_FDR"="grey", "n_sigFDR" = "red"),
                     labels = c("Permutations", "Permutations FDR", "RegEvol")) +
  scale_y_continuous(name = "Nb Directional", limits = c(0, max(df_plot$value))) +
  scale_x_continuous(breaks = seq(min(df_plot$cutoff), max(df_plot$cutoff), by=1)) +
  labs(x = "N substitution cutoff", color = "" ) +
  theme_minimal(base_size = 11) + theme( legend.position = c(leg_x,0.65))

df_plot <- reshape2::melt(sweep, id.vars = "cutoff", measure.vars = c("perc_perm", "perc_perm_FDR", "percFDR"))
df_plot$variable <- factor(df_plot$variable, levels = c("perc_perm", "perc_perm_FDR", "percFDR"))

p_sweep_prop <- ggplot(df_plot, aes(x = cutoff, y = value, color = variable, group = variable)) +
  geom_line(size = 1) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("perc_perm" = "black", "perc_perm_FDR"="grey", "percFDR" = "red"),
                     labels = c("Permutations", "Permutations FDR", "RegEvol")) +
  scale_y_continuous(name = "% Directional", limits = c(0, max(df_plot$value))) +
  scale_x_continuous(breaks = seq(min(df_plot$cutoff), max(df_plot$cutoff), by=1)) +
  labs(x = "N substitution cutoff", color = "" ) +
  theme_minimal(base_size = 11) + theme( legend.position ="none")

# Ascertainment Bias
# Conclusion FDR
maxLL$FDR_null_purif <- p.adjust(maxLL$p_value_null_purif, method="fdr")
maxLL$FDR_null_pos <- p.adjust(maxLL$p_value_null_pos, method="fdr")
maxLL$FDR_purif_pos <- p.adjust(maxLL$p_value_purif_pos, method="fdr")

maxLL$FDR_Conclusion <- "Neutral"
maxLL$FDR_Conclusion <- ifelse(maxLL$FDR_null_purif <= treshold, "Stabilizing", maxLL$FDR_Conclusion)
maxLL$FDR_Conclusion <- ifelse(maxLL$FDR_purif_pos <= treshold, ifelse(maxLL$AlphaPos>maxLL$BetaPos,"Directional (+)", "Directional (-)"), maxLL$FDR_Conclusion)
maxLL$FDR_Conclusion <- factor(maxLL$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
maxLL$FDR_binary <- ifelse(grepl("Directional", maxLL$FDR_Conclusion), 1, 0) 

# Conclusion qvalues
maxLL$qval_null_purif <- qvalue(maxLL$p_value_null_purif)$q
maxLL$qval_null_pos <- qvalue(maxLL$p_value_null_pos)$q
maxLL$qval_purif_pos <- qvalue(maxLL$p_value_purif_pos)$q

maxLL$qval_Conclusion <- "Neutral"
maxLL$qval_Conclusion <- ifelse(maxLL$qval_null_purif <= treshold, "Stabilizing", maxLL$qval_Conclusion)
maxLL$qval_Conclusion <- ifelse(maxLL$qval_purif_pos <= treshold, ifelse(maxLL$AlphaPos>maxLL$BetaPos,"Directional (+)", "Directional (-)"), maxLL$qval_Conclusion)
maxLL$qval_Conclusion <- factor(maxLL$qval_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
maxLL$qval_binary <- ifelse(grepl("Directional", maxLL$qval_Conclusion), 1, 0) 

par(mfrow=c(2,2))
p1_fdr <- ggplot(data.frame(p = maxLL$FDR_null_purif), aes(x = p)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "", x = "FDR", y = "Count") +
  theme_minimal(base_size = 11)

p2_fdr <- ggplot(data.frame(p = maxLL$FDR_purif_pos), aes(x = p)) +
  geom_histogram(bins = 50, fill = "salmon", color = "black") +
  labs(title = "", x = "FDR", y = "Count") +
  theme_minimal(base_size = 11)

p3_fdr <- ggplot(data.frame(p = maxLL$FDR_null_pos), aes(x = p)) +
  geom_histogram(bins = 50, fill = "lightgreen", color = "black") +
  labs(title = "", x = "FDR", y = "Count") +
  theme_minimal(base_size = 11)

p_sub <- ggplot(data.frame(p = maxLL$Nmut), aes(x = p)) +
  geom_histogram(bins = 50, fill = "grey", color = "black") +
  labs(title = "", x = "Substitutions", y = "Count") +
  theme_minimal(base_size = 11) + xlim(0, subMax)


all <- ((p1 | p2 | p3) /
    (p1_fdr | p2_fdr | p3_fdr) /
    (p_sub | p_sweep | p_sweep_prop)) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 14, face = "bold", hjust = -0.5), plot.tag.position = c(0, 1))
print(all)

ggsave(paste0(pathFigures, "/SupFig_pval_sub_effects_", sp, ".pdf"), all, width=9, height=9)
ggsave(paste0(pathFigures, "/SupFig_pval_sub_effects_", sp, ".png"), all, width=9, height=9)

################################################################################ 
# Ascertainment Bias
################################################################################
# Merge with rand
maxLL <- maxLL[which(maxLL$Nmut >= 4),]
rand <- read.table(file, h=T)
rownames(rand) <- rand$ID
rand <- rand[match(rownames(maxLL), rand$ID),]
rand$qval_binary <- maxLL$qval_binary
rand$FDR_binary <- maxLL$FDR_binary

rand <- rand[which(!is.na(rand$SVM)),]

rand$classSVM <- cut(rand$SVM, breaks=quantile(rand$SVM, probs = seq(0, 1, 0.1)), include.lowest = T)
rand$signif.bin <- ifelse(rand$pval.high <= 0.01, 1, 0)
rand$FDR <- p.adjust(rand$pval.high, method="fdr")
rand$signif.FDR.bin <- ifelse(rand$FDR <= treshold , 1, 0)

# Conclusion FDR corrected 
maxLL$FDR_null_purif <- p.adjust(maxLL$p_value_null_purif, method="fdr")
maxLL$FDR_null_pos <- p.adjust(maxLL$p_value_null_pos, method="fdr")
maxLL$FDR_purif_pos <- p.adjust(maxLL$p_value_purif_pos, method="fdr")

maxLL$FDR_Conclusion <- "Neutral"
maxLL$FDR_Conclusion <- ifelse(maxLL$FDR_null_purif <= treshold, "Stabilizing", maxLL$FDR_Conclusion)
maxLL$FDR_Conclusion <- ifelse(maxLL$FDR_purif_pos <= treshold, ifelse(maxLL$AlphaPos>maxLL$BetaPos,"Directional (+)", "Directional (-)"), maxLL$FDR_Conclusion)
maxLL$FDR_Conclusion <- factor(maxLL$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
maxLL$FDR_binary <- ifelse(grepl("Directional", maxLL$FDR_Conclusion), 1, 0) 

rand$FDR_binary <- maxLL$FDR_binary[match(rand$ID, rownames(maxLL))]
################################################################################
# Compute proportions for each fraction
fractions <- seq(0, 0.9, 0.1)
labels <- paste0(100 - fractions * 100, "%")

prop_perm <- c()
prop_perm_FDR <- c()
prop_maxll_FDR <- c()
for (q in fractions) {
  sub <- rand[which(rand$SVM >= quantile(rand$SVM, q)),]
  
  prop_perm <- c(prop_perm, sum(sub$signif.bin)/nrow(sub))
  prop_perm_FDR <- c(prop_perm_FDR, sum(sub$signif.FDR.bin)/nrow(sub))
  
  prop_maxll_FDR <- c(prop_maxll_FDR, sum(sub$FDR_binary, na.rm=T)/nrow(sub))
}

# Combine into a data.frame for ggplot
df_plot <- data.frame(
  Fraction = rep(labels, 3),
  Proportion = c(prop_perm, prop_perm_FDR, prop_maxll_FDR),
  Test = rep(c("Permutation", "Permutation_adjust", "RegEvol"), each = length(labels))
)

################################################################################
# Ensure proper x-axis order (descending fraction)
df_plot$Fraction <- factor(df_plot$Fraction, levels = labels)

p_bias <- ggplot(df_plot, aes(x = Fraction, y = Proportion, group = Test, color = Test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Permutation" = "black", "Permutation_adjust" = "darkgrey", "RegEvol" = "red"),
                     labels = c("Permutations", "Permutations FDR", "RegEvol")) +
  labs(title = "",  x = "Fraction of dataset (highest SVM)",
       y = "Prop. directional peaks", color = "" ) +
  ylim(0, max(df_plot$Proportion)*1.2) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "inside", legend.position.inside = c(0.2, 0.8))

  
plot(p_bias)
ggsave(paste0(pathFigures, "/SupFig_Ascertainment_Bias_", sp, ".pdf"), p_bias, width=6, height=6)
ggsave(paste0(pathFigures, "/SupFig_Ascertainment_Bias_", sp, ".png"), p_bias, width=6, height=6)
################################################################################ 

