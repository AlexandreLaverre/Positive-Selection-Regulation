library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

pathFigure = "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"

################################################################################      
# Get Data for plotting from Figure3.R script outputs
# Combine all simulations
cols = c("p_value_null_purif", "p_value_purif_pos", "p_value_null_pos",
         "FDR_null_purif", "FDR_purif_pos", "FDR_null_pos")

all_data <- bind_rows(
  all_pos %>% select(all_of(cols)) %>% mutate(simulation = "Directional Simulation"),
  all_null %>% select(all_of(cols)) %>% mutate(simulation = "Neutral Simulation"),
  all_stab %>% select(all_of(cols)) %>% mutate(simulation = "Stabilising Simulation"))

# pivot longer
all_tidy <- all_data %>%
  pivot_longer(
    cols = -simulation,
    names_to = c("metric", "comparison"),
    names_pattern = "(p_value|FDR)_(.*)",
    values_to = "value")

# Order and names
all_tidy$comparison <- factor(all_tidy$comparison, levels = c("null_purif", "purif_pos", "null_pos"), 
                              labels = c("Neutral vs Stabilising", "Stabilising vs Directional", "Neutral vs Directional"))
all_tidy <- all_tidy %>% mutate(metric = factor(metric, levels = c("p_value", "FDR")))

################################################################################      
# Plot 1: only p-values
plot_pval <- all_tidy %>%
  filter(metric == "p_value") %>%
  ggplot(aes(x = value, fill = comparison)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 color = "black", alpha = 0.7, position = "identity") +
  facet_wrap(~simulation, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Neutral vs Stabilising" = "skyblue",
    "Stabilising vs Directional" = "salmon",
    "Neutral vs Directional" = "lightgreen")) +
  labs(x = "p-values", y = "Density", fill = "LRT Comparisons") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    strip.background = element_rect(fill = "gray95", color = "black"),
    strip.text = element_text(face = "bold"))

# Plot 2: only FDR
plot_fdr <- all_tidy %>%
  filter(metric == "FDR") %>%
  ggplot(aes(x = value, fill = comparison)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 color = "black", alpha = 0.7, position = "identity") +
  facet_wrap(~simulation, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Neutral vs Stabilising" = "skyblue",
    "Stabilising vs Directional" = "salmon",
    "Neutral vs Directional" = "lightgreen")) +
  labs(x = "FDR", y = "Density", fill = "LRT Comparisons") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "black", size = 1),
    strip.background = element_rect(fill = "gray95", color = "black"),
    strip.text = element_text(face = "bold"))

combined <- plot_pval / plot_fdr + plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold"))
combined

# Save plots
ggsave(paste0(pathFigure, "/SupFig4_simul.png"), combined, width = 10, height = 8, dpi = 320)
ggsave(paste0(pathFigure, "/SupFig4_simul.pdf"), combined, width = 10, height = 8, dpi = 320)

################################################################################