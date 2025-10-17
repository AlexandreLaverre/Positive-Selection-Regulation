library(ggplot2)
library(ggExtra)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, "_filtered.Rds")
all_MLE_list <- readRDS(AllMerged)
all_MLE <- do.call(rbind, all_MLE_list)
FDR <- subset(all_MLE, all_MLE$FDR_Conclusion %in% c("Directional (+)", "Directional (-)"))
other <- subset(all_MLE, !(all_MLE$FDR_Conclusion %in% c("Directional (+)", "Directional (-)")))

# Other
Sub_exp_other = as.data.frame(tapply(other$Nmut, other$exp, function(x) mean(x)))
SNP_exp_other = as.data.frame(tapply(other$NbSNP, other$exp, function(x) mean(x)))
colnames(Sub_exp_other) = "Sub"
colnames(SNP_exp_other) = "SNP"
SNP_exp_other$Sub = Sub_exp_other$Sub 

# Directional
Sub_exp_FDR = as.data.frame(tapply(FDR$Nmut, FDR$exp, function(x) mean(x)))
SNP_exp_FDR = as.data.frame(tapply(FDR$NbSNP, FDR$exp, function(x) mean(x)))
colnames(Sub_exp_FDR) = "Sub"
colnames(SNP_exp_FDR) = "SNP"
SNP_exp_FDR$Sub = Sub_exp_FDR$Sub 

# Plot
SNP_exp_FDR$data <- "Directional"
SNP_exp_other$data <- "Other"
all <- rbind(SNP_exp_FDR, SNP_exp_other)

# Run Wilcoxon tests (results printed in terminal)
cat("Wilcoxon SNP (Directional < Other):\n")
print(wilcox.test(SNP ~ data, data = all, alternative = "less"))

cat("\nWilcoxon Sub (Directional > Other):\n")
print(wilcox.test(Sub ~ data, data = all, alternative = "greater"))

# Fit linear models
lm_dir <- lm(Sub ~ SNP, data = subset(all, data == "Directional"))
lm_other <- lm(Sub ~ SNP, data = subset(all, data == "Other"))

# Extract coefficients
coef_dir <- coef(lm_dir)
coef_other <- coef(lm_other)

# Format equations y = ax + b
eq_dir <- paste0("y = ", round(coef_dir[2], 3), "x + ", round(coef_dir[1], 2))
eq_other <- paste0("y = ", round(coef_other[2], 3), "x + ", round(coef_other[1], 2))

# Plot
p <- ggplot(all, aes(x = SNP, y = Sub, color = data)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(name = "", values = c("Directional" = "red", "Other" = "blue")) +
  labs(y = "Substitutions", x = "SNPs") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  annotate("text", x = max(all$SNP) * 0.6, y = min(all$Sub) + 7,
           label = eq_dir, color = "red", hjust = 0, size = 5) +
  annotate("text", x = max(all$SNP) * 0.6, y = min(all$Sub) + 5,
           label = eq_other, color = "blue", hjust = 0, size = 5)

# Add marginal distributions
p_marginal <- ggMarginal(p, type = "density", groupFill = TRUE, alpha = 0.4, margins = "both")
p_marginal

# Save plot
ggsave(filename = paste0("/Users/alaverre/Documents/Detecting_positive_selection/final_figures/Sup_Fig_SNP_vs_Sub_drosophila.png"),
       plot = p_marginal, width = 8, height = 6, dpi = 300)


