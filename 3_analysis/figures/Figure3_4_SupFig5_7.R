library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(ggtext)
par(xpd = TRUE)

sp="drosophila"
TF="Ni12/CTCF"
path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/", sp, "/", TF, "/")
pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"

bin = "exact_ranked_50"
selection = c("stabilising", "positive", "neutral")
maj = c("Stabilizing", "Positive", "Random")
names(maj) <- selection
maxSub = 150
maxLen = 1000
fdr_threshold = 0.01
pval_threshold = 0.01

params_pos=c(5, 10, 25, 50, 100) # pos
params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
params_purif=c("0.0", 0.1, 0.2, 0.5, "1.0", "2.0", "3.0", "5.0", "10.0")

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
    }else if (type =="stab"){
      data = paste0("beta_0.0Null_", param, "Stab_25Pos_")
    }
    
    obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
    deltas <- read.table(paste0(path, "/deltas/simul_", data, addExtreme, sel, "_observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
    row.names(deltas) <- deltas$seq_name
    all_delta[[paste0(sel, data)]] <- deltas
    
    simul <- read.csv(paste0(path, "/Tests/MLE_summary_simulated_", data, addExtreme, sel, "_", bin, "bins_threshold_0.01.csv"), h=T, row.names = 1)
    simul$Conclusion <- as.factor(simul$Conclusion)
    simul$Conclusion <- factor(simul$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
    simul$Conclusion_bin <- ifelse(grepl("Directional", simul$Conclusion), 1, 0)
    
    # FDR correction
    simul$FDR_null_pos <- p.adjust(simul$p_value_null_pos, method="fdr")
    simul$FDR_null_purif <- p.adjust(simul$p_value_null_purif, method="fdr")
    simul$FDR_purif_pos <- p.adjust(simul$p_value_purif_pos, method="fdr")
    
    simul$Conclusion_FDR <- "Neutral"
    simul$Conclusion_FDR <- ifelse(simul$FDR_null_purif <= fdr_threshold, "Stabilising", simul$Conclusion_FDR)
    simul$Conclusion_FDR <- ifelse(simul$FDR_purif_pos <= fdr_threshold, ifelse(simul$AlphaPos>simul$BetaPos,"Directional (+)", "Directional (-)"), simul$Conclusion_FDR)
    simul$Conclusion_FDR <- factor(simul$Conclusion_FDR, levels=c("Directional (+)", "Directional (-)", "Stabilising", "Neutral"))
    simul$Conclusion_FDR_bin <- ifelse(grepl("Directional", simul$Conclusion_FDR), 1, 0)
                                       
    simul$ID <- rownames(simul)
    simul$Scenario <- maj[sel]
    simul$param <- param
    simul$SVM <- deltas[rownames(simul), "SVM"]
    simul$deltaSVM <- deltas[row.names(simul),]$deltaSVM
    
    # Permut
    permut <- read.table(paste0(path, "/Tests/PosSelTest_deltaSVM_10000permutations_simulation_", data, addExtreme, sel, ".txt"), h=T, row.names = 1)
    permut <- permut[rownames(simul),]
    permut$signif_bin <- ifelse(permut$pval.two.tailed <= pval_threshold, 1, 0)
    permut$signif <- as.factor(ifelse(permut$pval.two.tailed > pval_threshold, "Neutral", ifelse(permut$deltaSVM >0, "Directional (+)", "Directional (-)")))
    permut$signif <- factor(permut$signif, levels=c("Directional (+)", "Directional (-)", "Neutral"))
    
    permut$FDR <- p.adjust(permut$pval.two.tailed, method="fdr")
    permut$signif_FDR_bin <- ifelse(permut$FDR <= fdr_threshold, 1, 0)
    permut$signif_FDR <- as.factor(ifelse(permut$FDR > fdr_threshold, "Neutral", ifelse(permut$deltaSVM >0, "Directional (+)", "Directional (-)")))
    permut$signif_FDR <- factor(permut$signif_FDR, levels=c("Directional (+)", "Directional (-)", "Neutral"))

    simul <- cbind(simul, permut[,c("pval.two.tailed", "signif", "signif_bin","signif_FDR", "signif_FDR_bin")])
    all_simul[[paste0(param)]] <- simul
  }
  
  all <- do.call(rbind, all_simul)
  
  return(all)
}

plotLines <- function(all, param, xlab="Directional Selection Strenght (alpha)", max=1, ylab="True Positive") {
  prop_TP_delta = prop.table(table(all$Conclusion_FDR, all[[param]]), margin=2)
  prop_TP_delta_signif = prop.table(table(all$signif, all[[param]]), margin=2)
  #prop_TP_delta_signif_FDR = prop.table(table(all$signif_FDR, all[[param]]), margin=2)
  param_order <- colnames(prop_TP_delta)

  # Convert to long format for ggplot
  if (all$Scenario[1] %in% c("Neutral", "Random")){ylab <- "False Positive Rate"}
    
  if (all$Scenario[1] == "Stabilizing" & ylab == "True Positive"){
    prop_total = c(prop_TP_delta[3,], rep(0, ncol(prop_TP_delta)))
  }else{
    prop_total = c(prop_TP_delta[1,] + prop_TP_delta[2,], 
                   prop_TP_delta_signif[1,] + prop_TP_delta_signif[2,])
  }

  df_plot <- data.frame(
    param = rep(colnames(prop_TP_delta), 2),
    prop_total = prop_total,
    Type = rep(c("RegEvol", "Permutations"), each = ncol(prop_TP_delta))
  )
  df_plot$Type <- factor(df_plot$Type, levels = c("RegEvol", "Permutations"))
  
  # Make param a factor to enforce order
  df_plot$var <- factor(df_plot$param, levels = param_order)
  df_plot$var <- as.numeric(df_plot$param)
  
  #if (param == "classdeltaSVM") {
  #  mid_labels <- sapply(levels(df_plot$var), function(x) {
  #    nums <- as.numeric(unlist(regmatches(x, gregexpr("-?\\d+\\.?\\d*", x))))
  #    round(mean(nums), 1)})
  #  
  #  df_plot$var <- factor(rep(as.character(mid_labels),2), levels = mid_labels)
  #}
  
  # Plot
  ggplot(df_plot, aes(x = var, y = prop_total, group = Type, color = Type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("RegEvol" = "red", "Permutations" = "black")) +
    labs(x = xlab, y = ylab, color = "") +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "inside",
          legend.position.inside = c(0.2, 0.8)) +
    ylim(0, max)  # same as in base R
}

################################################################################
# Figure 3 C-F
addExtreme = "0.0Purif_addExtreme0_" #0.0Purif_addExtreme0_
params_pos=c(11,12,13,14)
all_pos = get_all_simul(params_pos, type="pos", sel="positive")
breaks <- quantile(all_pos$deltaSVM, probs = seq(0, 1, 0.05), na.rm = TRUE)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
all_pos$classdeltaSVM <- cut(all_pos$deltaSVM, breaks=breaks, include.lowest = T, labels = midpoints)

p1 <- plotLines(all_pos, "Nmut", xlab="Substitution number", max=1)
p2 <- plotLines(all_pos, "classdeltaSVM", xlab=expression(Delta*"SVM"), max=1)

params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
all_prop_null = get_all_simul(params_null, type="null", sel="positive")
p3 <- plotLines(all_prop_null, "param", xlab="Proportion of Neutral Substitution", max=1)

params_pos=c(2, 3, 4, 5, 10, 25, 50, 75, 100)
all = get_all_simul(params_pos, type="pos", sel="positive")
p4 <- plotLines(all, "param", xlab=expression("Directional Selection Strength (" * alpha * ")"), max=1)


# Remove legend from all except the first
p1 <- p1 + theme(legend.position="inside", legend.position.inside = c(0.8, 0.3))
p2 <- p2 + theme(legend.position="none", axis.title.y = element_blank())
p3 <- p3 + theme(legend.position="none")
p4 <- p4 + theme(legend.position="none", axis.title.y = element_blank())

# Combine the plots side by side or in a grid
combined <- (p1 | p2) / (p3 | p4)   # 2x2 grid
combined

if (length(params_pos) == 1){
  sup5 <- p1 + p2 + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
  ggsave(paste0(pathFigure, "SuppFigure5_simul_positive.pdf"), plot=sup5, width = 9.5, height = 5, dpi = 320)
  ggsave(paste0(pathFigure, "SuppFigure5_simul_positive.png"), plot=sup5, width = 9.5, height = 5, dpi = 320)
}

################################################################################
# Stab
params_stab=c(11, 12, 13, 14)
params_stab=c(2, 3, 4, 5, 10, 25, 50, 75)
params_stab=c("40.0", 40.5, 41, 41.5, 42, 42.5, 43, 43.5, 44, 44.5, 45, 45.5, 46, 46.5, 47, 47.5, 48, 48.5, 49, 49.5)
#params_stab=c(40, 41, 42, 43, 44, 45,  46,  47,  48, 49) #41.5

all_stab = get_all_simul(params_stab, type="stab", sel="stabilising")
breaks <- quantile(all_stab$deltaSVM, probs = seq(0, 1, 0.01), na.rm = TRUE)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
all_stab$classdeltaSVM <- cut(all_stab$deltaSVM, breaks=breaks, labels = midpoints, include.lowest = T)

p_stab1 <- plotLines(all_stab, "Nmut", xlab="Substitution number", max=1)
p_stab_mut_false <- plotLines(all_stab, "Nmut", xlab="Nb Substitution (stabilising selection)", ylab="False Positive Rate", max=0.05) + theme(legend.position="none")
p_stab2 <- plotLines(all_stab, "classdeltaSVM", xlab=expression(Delta*"SVM"), max=1)
p_stab_svm_false <- plotLines(all_stab, "classdeltaSVM", xlab=expression(Delta*"SVM (stabilising selection)"), ylab="False Positive Rate", max=0.4)
p_stab_svm_false <- p_stab_svm_false + theme(legend.position="none", axis.title.y = element_blank())
  
params_stab=c(2, 3, 4, 5, 10, 25, 50, 75)
all_a_stab = get_all_simul(params_stab, type="stab", sel="stabilising")
p_stab4 <- plotLines(all_a_stab, "param", xlab=expression("Stabilising Selection Strength ("*alpha*"="*beta*")"), max=1)

params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
all_prop_stab = get_all_simul(params_null, type="null", sel="stabilising")
p_stab3 <- plotLines(all_prop_stab, "param", xlab="Proportion of Neutral Substitution", max=1)

p_stab1 <- p_stab1 + theme(legend.position="inside", legend.position.inside = c(0.8, 0.3))
p_stab2 <- p_stab2 + theme(legend.position="none", axis.title.y = element_blank())
p_stab3 <- p_stab3 + theme(legend.position="none")
p_stab4 <- p_stab4 + theme(legend.position="none", axis.title.y = element_blank())
combined_stab <- (p_stab1 | p_stab2) / (p_stab3 | p_stab4) +
                  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
combined_stab
ggsave(paste0(pathFigure, "SupFigure_simul_stabilising.pdf"), width = 9.5, height = 7, dpi = 320)
ggsave(paste0(pathFigure, "SupFigure_simul_stabilising.png"), width = 9.5, height = 7, dpi = 320)


################################################################################
###### Stacked barplots
params_null=c(11, 12, 13, 14)
addExtreme = "0.0Purif_addExtreme0_" #0.0Purif_addExtreme0_
all_null = get_all_simul(params_null, type="pos", sel="neutral")
all_null$classdeltaSVM <- cut(all_null$deltaSVM, breaks=quantile(all_null$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)

all_pos$Scenario  <- ifelse(all_pos$AlphaPos>all_pos$BetaPos, "Directional(+)", "Directional(-)")
all_null$Scenario <- "Neutral"
all_stab$Scenario <- "Stabilising"

# Rename for clarity
df <- rbind(
  all_stab %>% mutate(Method="RegEvol", InferredModel=Conclusion_FDR),
  all_pos  %>% mutate(Method="RegEvol", InferredModel=Conclusion_FDR),
  all_null %>% mutate(Method="RegEvol", InferredModel=Conclusion_FDR),
  all_stab %>% mutate(Method="Permutations", InferredModel=signif),
  all_pos  %>% mutate(Method="Permutations", InferredModel=signif),
  all_null %>% mutate(Method="Permutations", InferredModel=signif)
)

# Make sure both are factors with consistent levels
levels_order <- c("Directional (+)", "Directional (-)", "Neutral", "Stabilising")
df$InferredModel <- factor(df$InferredModel, levels = levels_order)

# Reshape long for plotting
plot_df <- df %>%
  group_by(Method, Scenario, InferredModel) %>%
  summarise(N = n(), .groups = "drop_last") %>%
  mutate(Proportion = N / sum(N))

plot_df$Method <- factor(plot_df$Method, levels = c("RegEvol", "Permutations"))
cb_palette <- c(
  "Directional (+)" = "orange",
  "Directional (-)" = "navy",
  "Stabilising"     = "#009E73",
  "Neutral"         = "lightgrey"
)

# Plot
pbars <- ggplot(plot_df, aes(x = Scenario, y = Proportion, fill = InferredModel)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  facet_wrap(~Method) +
  scale_fill_manual(values = cb_palette) +
  labs(x = "Simulated evolutionary scenario", y = "Proportion", fill = "Inferred model") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.5, 1.3),
    legend.direction = "horizontal",
    plot.margin = margin(t = 30, r = 5, b = 5, l = 5),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )


# Combine rows
final_plot <- (pbars) / 
  (p1 | p2) / 
  (p3 | p4) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))

# Display
final_plot

ggsave(paste0(pathFigure, "Figure3.pdf"), final_plot, width = 11, height = 11)
ggsave(paste0(pathFigure, "Figure3.png"), final_plot, width = 11, height = 11)

################################################################################
# Figure 4 
# False positive on delta SVM
signed_log <- function(x, base = 10) {
  sign(x) * log1p(abs(x)) / log(base)
}
all_null$deltaSVM_log <- signed_log(all_null$deltaSVM)
breaks <- quantile(all_null$deltaSVM_log, probs = seq(0, 1, 0.01), na.rm = TRUE)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
all_null$classdeltaSVM <- cut(all_null$deltaSVM_log, breaks=breaks, include.lowest = T, labels=midpoints)

p_null_svm_false <- plotLines(all_null, "classdeltaSVM", xlab=expression(Delta*"SVM (random drift)"), max=0.4, ylab="False Positive Rate")
p_null_svm_false <- p_null_svm_false + theme(legend.position="inside", legend.position.inside = c(0.3, 0.8))
p_null_svm_false
# Ascertainment Bias
all_null$new_SVM <- all_null$SVM + all_null$deltaSVM

# Get subset of data with top X% SVM
fractions <- seq(0, 0.9, 0.1)
labels <- paste0(100 - fractions * 100, "%")
prop_perm <- c()
prop_maxll_FDR <- c()
for (q in fractions) {
  sub <- all_null[which(all_null$new_SVM >= quantile(all_null$new_SVM, q)),]
  prop_perm <- c(prop_perm, sum(sub$signif_bin)/nrow(sub))
  prop_maxll_FDR <- c(prop_maxll_FDR, sum(sub$Conclusion_FDR_bin, na.rm=T)/nrow(sub))
}

# Combine into a data.frame for ggplot
df_plot <- data.frame(
  Fraction = rep(labels, 2),
  Proportion_SVM = c(prop_perm,  prop_maxll_FDR),
  Test = rep(c("Permutation", "RegEvol"), each = length(labels))
)

df_plot$Fraction <- factor(df_plot$Fraction, levels = labels)

p_bias <- ggplot(df_plot, aes(x = Fraction, y = Proportion_SVM, group = Test, color = Test)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Permutation" = "black", "RegEvol" = "red")) +
  labs( title = "", x = "Fraction of dataset (highest SVM)", y = "False Positive Rate", color = "" ) +
  ylim(0, 0.05) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", axis.title.y = element_blank())

fig4 <- p_null_svm_false + p_stab_svm_false + p_stab_mut_false + p_bias  + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
fig4
ggsave(paste0(pathFigure, "Figure4.pdf"), fig4, width = 9.5, height = 7, dpi = 320)
ggsave(paste0(pathFigure, "Figure4.png"), fig4, width = 9.5, height = 7, dpi = 320)

