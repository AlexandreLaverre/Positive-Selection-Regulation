library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
par(xpd = TRUE)

sp="drosophila"
TF="Ni12/CTCF"
path = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/", sp, "/", TF, "/")
pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

bin = "exact_ranked_50"
selection = c("stabilising", "positive", "neutral")
maj = c("Stabilizing", "Positive", "Random")
names(maj) <- selection
maxSub = 150
maxLen = 1000

params_pos=c(2, 5, 10, 25, 50, 100) # pos
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

plotLines <- function(all, param, xlab="Directional Selection Strenght (alpha)", ylab="Proportion True Positive", max=1) {
  prop_TP_delta = prop.table(table(all$Conclusion, all[[param]]), margin=2)
  prop_TP_delta_signif = prop.table(table(all$signif, all[[param]]), margin=2)
  param_order <- colnames(prop_TP_delta)

  # Convert to long format for ggplot
  df_plot <- data.frame(
    param = rep(colnames(prop_TP_delta), 2),
    prop_total = c(prop_TP_delta[1,] + prop_TP_delta[2,], 
                   prop_TP_delta_signif[1,] + prop_TP_delta_signif[2,]),
    Type = rep(c("RegEvol", "Permutations"), each = ncol(prop_TP_delta))
  )
  df_plot$Type <- factor(df_plot$Type, levels = c("RegEvol", "Permutations"))
  
  # Make param a factor to enforce order
  df_plot$var <- factor(df_plot$param, levels = param_order)
  
  if (param == "classdeltaSVM") {
    mid_labels <- sapply(levels(df_plot$var), function(x) {
      nums <- as.numeric(unlist(regmatches(x, gregexpr("-?\\d+\\.?\\d*", x))))
      round(mean(nums), 1)})
    
    df_plot$var <- factor(rep(as.character(mid_labels),2), levels = mid_labels)
  }
  
  # Plot
  ggplot(df_plot, aes(x = var, y = prop_total, group = Type, color = Type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("RegEvol" = "black", "Permutations" = "red")) +
    labs(x = xlab, y = ylab, color = "") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.2)) +
    ylim(0, max)  # same as in base R
}

################################################################################
# Figure 3 C-F
addExtreme = "0.0Purif_addExtreme0_"
params_pos=c(25, 50, 100)
all = get_all_simul(params_pos, type="pos", sel="positive")
all$classdeltaSVM <- cut(all$deltaSVM, breaks=quantile(all$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)

p1 <- plotLines(all, "Nmut", xlab="Substitution number", max=1)
p2 <- plotLines(all, "classdeltaSVM", xlab=expression(Delta*"SVM class"), max=1) +
      scale_x_discrete(labels = function(x) {
      labs <- as.character(x)
      labs[seq_along(labs) %% 2 != 1] <- ""  # show every 2th label
      labs
    })

params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
all = get_all_simul(params_null, type="null", sel="positive")
p3 <- plotLines(all, "param", xlab="Proportion of Neutral Substitution", max=1)

params_pos=c(2, 5, 10, 25, 50, 100) 
all = get_all_simul(params_pos, type="pos", sel="positive")
p4 <- plotLines(all, "param", xlab=expression("Directional Selection Strength (" * alpha * ")"), max=1)


# Remove legend from all except the first
p1 <- p1
p2 <- p2 + theme(legend.position="none", axis.title.y = element_blank())
p3 <- p3 + theme(legend.position="none")
p4 <- p4 + theme(legend.position="none", axis.title.y = element_blank())

# Combine the plots side by side or in a grid
combined <- (p1 | p2) / (p3 | p4)   # 2x2 grid

# Display
combined

################################################################################
# Figure 4 A-D
params_null=c("0.0", 0.01, 0.05, 0.1, 0.25, 0.5)
addExtreme = "0.0Purif_addExtreme0_"
all = get_all_simul(params_null, type="null", sel="neutral")
all$classdeltaSVM <- cut(all$deltaSVM, breaks=quantile(all$deltaSVM, probs = seq(0, 1, 0.05)), include.lowest = T)

p1 <- plotLines(all, "classdeltaSVM", xlab=expression(Delta*"SVM class"),
          ylab="Proportion False Positive", max=0.4) +
          scale_x_discrete(labels = function(x) {labs <- as.character(x) 
          labs[seq_along(labs) %% 2 != 1] <- ""  # show every 2th label
          labs})


addExtreme = "0.0Purif_addExtreme1_"
params_null=c("0.0")
all = get_all_simul(params_null, type="null", sel="neutral")
p2 <- plotLines(all, "Nmut", xlab="Substitution number",
                ylab="Proportion False Positive", max=0.4) + theme(legend.position="none", axis.title.y = element_blank())

p1 + p2



addExtreme = ""
params_pos=c(2, 5, 10, 25, 50, 100) 
all = get_all_simul(params_pos, type="pos", sel="positive")
all$classSVM <- cut(all$SVM, breaks=quantile(all$SVM, probs = seq(0, 1, 0.1)), include.lowest = T)
plotLines(all, "classSVM", xlab="SVM class", ylab="Proportion True Positive", max=1) + theme(legend.position="none")


# Stab
params_null=c("0.0", 0.1, 0.25)
all = get_all_simul(params_null, type="null", sel="stabilising")
plotLines(all, "Nmut", xlab="Substitution number", max=1)
