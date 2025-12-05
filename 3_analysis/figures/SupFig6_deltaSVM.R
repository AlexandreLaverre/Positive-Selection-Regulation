library(data.table)
library(ggplot2)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"

species <- c("drosophila", "mouse", "human")  
TFs <- c("CTCF")                              
samples <- c("Ni12", "Schmidt12")

################################################################################
all_deltas_summary <- list()
for (sp in species) {
  for (samp in samples) {
    for (tf in TFs) {
      print(paste("Processing:", sp, samp, tf))
      file_all_deltas <- file.path(path, "NarrowPeaks", sp, samp, tf, "deltas", "ancestral_all_possible_deltaSVM.txt")
      
      if (!file.exists(file_all_deltas)) {next}
      deltas_all <- fread(file_all_deltas, header = TRUE, fill = TRUE)
      rownames(deltas_all) <- deltas_all[[1]]
      deltas_all[[1]] <- NULL
      
      vals <- as.numeric(as.matrix(deltas_all))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) {next}
      
      all_deltas_summary[[paste0(sp, "_", TF)]] <- data.table(species = sp, TF = tf, deltaSVM = vals)
      }}}

# Combine all into one data.table
all_deltas_summary <- rbindlist(all_deltas_summary, fill = TRUE)
median_df <- all_deltas_summary %>%
  group_by(species, TF) %>%
  summarise(median_val = median(deltaSVM, na.rm = TRUE), .groups = "drop")
################################################################################
# Plot histogram for all deltaSVMs
hist <- ggplot(all_deltas_summary, aes(x = deltaSVM)) +
        geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 0.2, 
                       fill = "gray90", color = "black") +
        geom_vline(xintercept = 0, color = "black", size = 0.8) +
        geom_vline(data = median_df, aes(xintercept = median_val), color = "red", size = 1) +
        geom_text(data = median_df, aes(x = median_val-3, y = 0.01, label = paste0("med=", round(median_val,2))),
                  color = "red", vjust = -0.5, size=5) +
        coord_cartesian(xlim = c(-7.5, 7.5)) +
        facet_grid(species ~ TF) +
        labs(x = expression(Delta*"SVM"), y = "Proportion") +
        theme_minimal(base_size = 16)

  
# Save the plot
ggsave(file.path(pathFigures, "SupFig6_deltaSVM.png"), width = 8, height = 8, dpi = 320)
ggsave(file.path(pathFigures, "SupFig6_deltaSVM.pdf"), width = 8, height = 8, dpi = 320)

################################################################################