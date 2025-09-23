library(ggplot2)
library(dplyr)
library(ggpubr)


path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"

exp_info <- read.table("/Users/alaverre/Documents/Detecting_positive_selection/results/peaks_calling/peaks_numbers.txt", h=T)
rownames(exp_info) <- paste(exp_info$species, exp_info$TF, sep = " ")

models <- readRDS(paste0(path, "all_models.rds"))

cormat <- round(cor(models, ), 2)

# Reorder the matrix: use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
melted_cormat <- melt(cormat, na.rm = TRUE)

tf_colors <- setNames(ggpubr::get_palette("npg", length(TFs)), TFs)


# Get unique TFs and assign colors
TFs <- unique(exp_info$TF)
tf_colors <- setNames(ggpubr::get_palette("npg", length(TFs)), TFs)

SampleName <- rownames(exp_info)
# Summarize melted_cormat to get one label per sample
# We'll keep only one label per sample (species) for axis
melted_cormat$species <- exp_info$species[match(melted_cormat$Var1, SampleName)]
melted_cormat$TF <- exp_info$TF[match(melted_cormat$Var1, SampleName)]

# Reorder factor to match the correlation matrix order
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = unique(melted_cormat$Var2))
melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = unique(melted_cormat$Var1))

# Create the heatmap
ggplot(melted_cormat, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name = "Pearson\nCorrelation") +
  theme_minimal() +
  # Only species labels, colored by TF
  scale_x_discrete(labels = melted_cormat$species) +
  scale_y_discrete(labels = melted_cormat$species) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 14, color = tf_colors[melted_cormat$TF]),
    axis.text.y = element_text(size = 14, color = tf_colors[melted_cormat$TF]),
    axis.title = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12) 
  ) + 
  geom_point(aes(x = +Inf, y = +Inf, color = TF), size = 10) +
  scale_color_manual(name = "TF", values = tf_colors)

