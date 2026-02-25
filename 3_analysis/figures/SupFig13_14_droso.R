# Analyse peaks droso
library(stringr)
library(data.table)
library(ROCR)
library(ggplot2)
library(dplyr)
library(patchwork)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
method="exact_ranked_ancestral"
fdr_threshold = 0.05

AllMerged = paste0(path, "allMLE_drosophila_", method, "_2sub_all_exp.Rds")
all_MLE_list <- readRDS(AllMerged)
all_MLE <- do.call(rbind, all_MLE_list)

make_hist <- function(data, value, xlab) {
  med <- median(data[[value]], na.rm = TRUE)
  ggplot(data, aes_string(x = value)) +
    geom_histogram(fill = "gray70", color = "black", bins = 50) +
    geom_vline(xintercept = med, color = "red", size = 1) +
    labs(x = xlab, y = "Count") +
    theme_minimal(base_size = 14)
}

pA <- make_hist(all_MLE[!duplicated(all_MLE$exp),], "original_sample_size", "Peaks per experiment")
pB <- make_hist(all_MLE, "original_length", "Peak length (bp)")
pC <- make_hist(all_MLE, "Nmut", "Substitution number")
pD <- make_hist(all_MLE, "deltaSVM", expression(Delta*SVM))


all <- ((pA | pB ) / (pC | pD)) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 14, face = "bold", hjust = -0.5), plot.tag.position = c(0, 1))
print(all)

ggsave(paste0(pathFigures, "/SupFig_drosophila_stats.pdf"), all, width=9, height=9, dpi=320)
ggsave(paste0(pathFigures, "/SupFig_drosophila_stats.png"), all, width=9, height=9, dpi=320)

################################################################################
# Developmental stage annotation
all_MLE$DevStage <- gsub("_.*", "", all_MLE$DevStage)
all_MLE[grepl("instarlarva", all_MLE$DevStage),]$DevStage <- "larva"
all_MLE$DevStage <- factor(all_MLE$DevStage, levels=c("embryonic", "larva", "prepupa", "pupa", "adult", "unknown"))
all_MLE <- all_MLE[which(all_MLE$original_sample_size > 1500 & all_MLE$AUC >0.8),]

# number of peaks per DevStage
peaks_per_stage <- all_MLE %>% group_by(DevStage) %>%  summarise(N_peaks = n())

# number of experiments per DevStage
exp_per_stage <- all_MLE %>% distinct(exp, DevStage) %>% group_by(DevStage) %>% summarise(N_exp = n())

df_counts <- peaks_per_stage %>% left_join(exp_per_stage, by = "DevStage")

p1 <- ggplot(df_counts, aes(x = DevStage, y = N_peaks)) +
  geom_col(fill = "steelblue", color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  geom_text(aes(label = paste0("Nexp=", N_exp)), vjust = -0.5, size = 4) +
  labs(x = "Developmental stage", y = "Number of peaks") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 14))

################################################################################
# Proportion of peaks under directional selection per DevStage
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_threshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_threshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

df_prop <- all_MLE %>% group_by(DevStage) %>%
  summarise(N_dir = sum(FDR_Conclusion %in% c("Directional (+)", "Directional (-)")), 
            N_peaks = n(), proportion = N_dir / N_peaks)

p2 <- ggplot(df_prop, aes(x = DevStage, y = proportion)) +
  geom_col(fill = "salmon", color = "black") +
  geom_text(aes(label = paste0("N=", N_dir)), vjust = -0.5, size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Developmental stage",  y = "Proportion Directional") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 14))

all <- (p1 / p2 ) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 14, face = "bold", hjust = -0.5), plot.tag.position = c(0, 1))
print(all)


ggsave(paste0(pathFigures, "/SupFig_drosophila_stages.pdf"), all, width=7, height=7, dpi=320)
ggsave(paste0(pathFigures, "/SupFig_drosophila_stages.png"), all, width=7, height=7, dpi=320)
