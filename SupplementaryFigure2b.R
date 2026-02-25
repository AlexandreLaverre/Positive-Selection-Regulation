library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(ggseqlogo)
library(patchwork)
library(tidyr)
library(cowplot)
library(readr)
library(progress)

path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
all_data <- readRDS("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/Correlation_SVM_Coverage_Quality.Rds")

sp="human"
sample="Wilson"
TFS=c("CEBPA")
tf_palette <- ggpubr::get_palette("npg", 6)

# Combine into one dataframe
combined <- bind_rows(lapply(names(all_data), function(data){
  df <- all_data[[data]]
  # Ensure no NAs
  df$class_coverage <- cut(df$Coverage, breaks=quantile(df$Coverage, probs=seq(0,1, by=0.1), na.rm=TRUE), labels=seq(1,10), include.lowest = TRUE)
  df$class_SignalStrengthScore <- cut(df$SignalStrengthScore, breaks=quantile(df$SignalStrengthScore, probs=seq(0,1, by=0.1), na.rm=TRUE), labels=seq(1,10), include.lowest = TRUE)
  df %>% mutate(dataset = data, sp = str_split_i(data, "_", 1), TF = str_split_i(data, "_", 2)) %>% filter(!is.na(cor_svm_P), !is.na(Coverage), !is.na(SignalStrengthScore))
}))


# Signal Strength Score
summary_signal <- combined %>%
  group_by(dataset, class_SignalStrengthScore) %>%
  summarise(median = median(cor_svm_P, na.rm = TRUE),
            se_median = 1.58 * IQR(cor_svm_P) / sqrt(n()), .groups = "drop") # se_median = notch from boxplot

# Coverage Score
summary_cov <- combined %>%
  group_by(dataset, class_coverage) %>%
  summarise(median = median(cor_svm_P, na.rm = TRUE),
            se_median = 1.58 * IQR(cor_svm_P) / sqrt(n()), .groups = "drop") 

# Reorder legend
summary_signal$dataset <- factor(summary_signal$dataset, levels = c("human_CEBPA", "mouse_HNF4A", "mouse_HNF6",  "human_FOXA1", "human_HNF4A", "drosophila_CTCF"))
summary_cov$dataset <- factor(summary_cov$dataset, levels = c("human_CEBPA", "mouse_HNF4A", "mouse_HNF6",  "human_FOXA1", "human_HNF4A", "drosophila_CTCF"))

p_signal <- ggplot(summary_signal, aes(x = class_SignalStrengthScore, y = median,  group = dataset, color = dataset)) +
  geom_line() + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = median - se_median, ymax = median + se_median),  width = 0.2, size = 0.6) +
  labs(x = "Peak Signal Strength (quantile)", y = "Base-level SVM–Coverage correlation", color = "Dataset") +
  scale_color_manual(values = tf_palette, name = "") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.text = element_text(size = 14), 
        legend.position = "bottom", legend.margin = margin(t = 5, r = 5, b = 5, l = 300)) +
  guides(color = guide_legend(nrow = 2))

p_coverage <- ggplot(summary_cov, aes(x = class_coverage, y = median,  group = dataset, color = dataset)) +
  geom_line() + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = median - se_median, ymax = median + se_median),  width = 0.2, size = 0.6) +
  labs(x = "Peak Total Coverage (quantile)", y = "Base-level SVM–Coverage correlation", color = "Dataset") +
  scale_color_manual(values = tf_palette, name = "Dataset") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none")


# Coefficients distribution 
human_CEBPA <- all_data[[paste0("human_CEBPA")]]
df_correlation <- data.frame(correlation = human_CEBPA$cor_svm_P, pval = human_CEBPA$cor_svm_P_pval)
df_correlation <- df_correlation[complete.cases(df_correlation), ]
df_correlation$significant <- df_correlation$pval < 0.01

C <- ggplot(df_correlation, aes(x = correlation, fill = significant)) +
  geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
  labs(x = "Base-level SVM–Coverage correlation", y = "Count", fill = "p < 0.01") +
  theme_minimal(base_size = 14) +
  theme(legend.position = c(0.2, 0.8))

# Combine all correlations into one dataframe for violin plot
human_CEBPA <- all_data[[paste0("human_CEBPA")]]
human_CEBPA <- human_CEBPA[complete.cases(human_CEBPA), ]
N = nrow(human_CEBPA)
data_full <- data.frame(name=rep(c("SVM", "ΔSVM", "PhastCons", "PhyloP"), each = N), group="All peaks",
                        value=c(human_CEBPA$cor_svm_P, human_CEBPA$cor_delta_P, human_CEBPA$cor_phast_P, human_CEBPA$cor_phylo_P))

human_CEBPA_filtered <- human_CEBPA[which(as.numeric(human_CEBPA$class_SignalStrengthScore) > 5 & as.numeric(human_CEBPA$class_coverage) > 5),]
N_filtered = nrow(human_CEBPA_filtered)

data_filtered <- data.frame(name=rep(c("SVM", "ΔSVM", "PhastCons", "PhyloP"), each = N_filtered), group="High quality peaks",
                   value=c(human_CEBPA_filtered$cor_svm_P, human_CEBPA_filtered$cor_delta_P, human_CEBPA_filtered$cor_phast_P, human_CEBPA_filtered$cor_phylo_P))

data_all <- rbind(data_full, data_filtered)

pearson <- ggplot(data_all, aes(name, value)) +
  geom_violin(aes(fill=name, alpha=group), position=position_dodge(0.9), width=1.1) +
  geom_boxplot(aes(alpha=group), position=position_dodge(0.9), width=0.12, color="black") +
  scale_alpha_manual(values=c("All peaks"=0.2,"High signal"=0.9)) +
  scale_fill_manual(values=c("SVM"="#7CAE00","ΔSVM"="orange","PhastCons"="navy","PhyloP"="red"), guide="none") +
  scale_y_continuous(breaks=seq(-1,1,0.25)) +
  scale_x_discrete(limits=c("SVM","ΔSVM","PhastCons","PhyloP")) +
  theme_minimal(base_size=14) +
  theme(legend.position="none", axis.text=element_text(size=14, colour="black")) +
  xlab("") + ylab("Base-level SVM–Coverage correlation")

pearson

final_plot <- (p_signal + p_coverage) / (C + pearson) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))

final_plot

ggsave(paste0(pathFigures, "SupFig2_correlation_peak_quality.png"), plot = final_plot, width = 10, height = 10, dpi = 320)
ggsave(paste0(pathFigures, "SupFig2__correlation_peak_quality.pdf"), plot = final_plot, width = 10, height = 10, dpi = 320)

