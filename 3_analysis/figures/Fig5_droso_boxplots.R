library(ggplot2)
library(dplyr)
library(patchwork)

pathFigure <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
################################################################################
# Merge and bin
all_genes <- merge(all_genes_original, evol_data, by = "Gene.name") %>%
  mutate(Freq_bin = cut(Freq, 
                        breaks = c(0, 10, 25, 100, 500, 1000, 5500), 
                        include.lowest = TRUE,
                        labels = c("0-10", "11-25", "26-100", "101-500", "501-1000", "1001-5500")))
# Spearman correlation
cor <- cor.test(all_genes$Freq, all_genes[[omega]], method = "spearman")

ylab = ifelse(omega == "meanlog2_mel", "Mean gene expression log2(TPM)", "log2 Fold Change") # 
# Plot: violin (soft, background) + boxplot (emphasized) + jitter (faint points)
p1 <- ggplot(all_genes, aes(x = Freq_bin, y = !!sym(omega))) +
  geom_violin(fill = "gray70", color = NA, alpha = 0.6, trim = TRUE) +       # violin soft, no border
  geom_jitter(width = 0.2, alpha = 0.03, size = 0.5, color = "black") +      # dots faint
  geom_boxplot(outlier.shape = NA, notch = TRUE,
               fill = "white", alpha = 0.7, width = 0.5,                     # box slightly transparent
               color = "black", size = 0.6) +                               # thinner border
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Nb peaks",
    y = ylab,
  ) +
  annotate("text", x=3.5, y = max(all_genes[[omega]], na.rm = TRUE) * 0.95,
           label = paste0("Spearman's rho = ", round(cor$estimate, 2),
                          "\n p-val = ", format.pval(cor$p.value, digits = 2)),
           size = 5, hjust = 0.5)


################################################################################
# Ratio Positive
# Prepare data
pos_genes$all <- all_genes[match(pos_genes$Gene.name, all_genes$Gene.name),]$Freq
pos_genes_filtered <- pos_genes %>% 
  filter(all >= 2) %>%
  mutate(ratio = Freq / all,
         ratio_bin = cut(ratio, breaks = c(0, 0.05, 0.1, 0.2, 1), include.lowest = TRUE))

# Spearman correlation
cor_res <- cor.test(pos_genes_filtered$ratio, pos_genes_filtered[[omega]], method = "spearman")

# Boxplot with points
p2 <- ggplot(pos_genes_filtered, aes(x = ratio_bin, y = !!sym(omega))) +
  geom_violin(fill = "gray70", color = NA, alpha = 0.6, trim = TRUE) +       # violin soft, no border
  geom_jitter(width = 0.1, alpha = 0.03, size = 0.5, color = "black") +   # points behind
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = "white", color = "black", alpha=0.7) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Ratio Positive", y = "") +
  annotate("text", x = 2.5, y = max(pos_genes_filtered[[omega]], na.rm = TRUE) * 0.95,
    label = paste0("Spearman's rho = ", round(cor_res$estimate, 2), 
                   "\n p-val = ", format.pval(cor_res$p.value, digits=2)),
    size = 5,  hjust = 0.5)

# Combine plots
final_plot <- p1 + p2 
print(final_plot)
# Save plot
ggsave(paste0(pathFigure, "SupFig_droso_", omega, "_", anat, ".png"), 
       plot = final_plot, width = 14, height = 7, dpi = 300)

################################################################################