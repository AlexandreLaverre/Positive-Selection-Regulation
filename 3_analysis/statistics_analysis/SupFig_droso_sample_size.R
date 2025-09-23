library(ggplot2)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
sp = "drosophila"

method="exact_ranked_ancestral" 
all_MLE_list <- readRDS(paste0(path, "allMLE_drosophila_", method, ".Rds"))
all_MLE <- do.call(rbind, all_MLE_list)
all_MLE$peakID <- unlist(lapply(all_MLE_list, rownames))
rownames(all_MLE) <- all_MLE$peakID

all_MLE$DevStage <- gsub("_.*", "", all_MLE$DevStage)
all_MLE[grepl("instarlarva", all_MLE$DevStage),]$DevStage <- "larva"
all_MLE$DevStage <- factor(all_MLE$DevStage, levels=c("embryonic", "larva", "prepupa", "pupa", "adult", "unknown"))


SampleSize <- tapply(all_MLE$SampleSize, as.factor(all_MLE$exp), mean)
df <- data.frame(exp = names(SampleSize), sample_size = as.numeric(SampleSize))
df$stage <- all_MLE[match(df$exp, all_MLE$exp), "DevStage"]
df <- df[order(df$stage), ]

# Create consistent color palette
palette <- rainbow(length(levels(df$stage)))
names(palette) <- levels(df$stage)

ggplot(df, aes(x = stage, y = sample_size, fill = stage)) +
  geom_jitter(shape = 21, color = "black", width = 0.25, size = 2, alpha=0.5) +
  geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  theme(
    legend.position = "none",              # remove legend
    axis.text = element_text(size = 14),   # axis tick text
    axis.title = element_text(size = 15)   # axis labels
  ) +
  labs(title="N=674 experiments", y = "Number of Peaks", x = "")
ggsave(paste0(pathFigures, "/figS1_droso_sample_size.pdf", width = 6, height = 4, useDingbats = FALSE))
