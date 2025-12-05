library("ROCR")
library("data.table")
library("progress")
library("ggplot2")
library("ggExtra")

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathFigure <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
exp.folders <- list.dirs(path, recursive = FALSE)

# All results
method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, "_5sub_all_exp.Rds")
all_MLE_list <- readRDS(AllMerged)

# Get Model performances
all_AUC <- data.frame(exp = rep(NA,length(exp.folders)), AUC = NA, SampleSize=NA)
pb <- progress_bar$new(format = "  Processing [:bar] :percent (:elapsed s) | ETA :eta", total = length(exp.folders), clear = FALSE, width = 60)
for (i in 1:length(exp.folders)){
  pb$tick()
  exp.path = exp.folders[i]
  exp <- basename(exp.path)
  if (exp == "log"){next}
  
  cv<-fread(paste0(exp.path, "/Model/", exp, ".cvpred.txt"))
  colnames(cv)<-c("position","prediction","real_state","cv_number")
  pred <- prediction(cv$prediction, cv$real_state) 
  auc_result <- performance(pred, measure = "auc")
  AUC = signif(unlist(slot(auc_result, "y.values")), digits=3)
  
  all_AUC$exp[i] <- exp
  all_AUC$AUC[i] <- AUC
  all_AUC$SampleSize[i] <- unique(all_MLE_list[[exp]]$original_sample_size)
  all_MLE_list[[exp]]$AUC <- AUC
}
# Save results
saveRDS(all_MLE_list, AllMerged)


################################################################################
# Plot
all_AUC <- all_AUC[!is.na(all_AUC$AUC),]
rownames(all_AUC) <- all_AUC$exp

# Base scatterplot
p_auc <- ggplot(all_AUC, aes(x = SampleSize, y = AUC)) +
  geom_point(alpha = 0.6, size = 2, color = "black") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed", linewidth = 0.8) +
  theme_minimal(base_size = 16) +
  labs(x = "Sample size", y = "AUC") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Add marginal histograms
p_auc_marginal <- ggMarginal(p_auc, type = "histogram", fill = "gray60", color = "black", bins = 100)
p_auc_marginal

# Save plot
ggsave(paste0(pathFigure, "/FigSup16_AUC_vs_SampleSize_drosophila.png"), plot = p_auc_marginal, width = 8, height = 6, dpi = 320)
ggsave(paste0(pathFigure, "/FigSup16_AUC_vs_SampleSize_drosophila.pdf"), plot = p_auc_marginal, width = 8, height = 6, dpi = 320)
################################################################################