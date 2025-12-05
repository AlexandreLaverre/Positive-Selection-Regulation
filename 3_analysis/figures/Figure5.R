# Analyse peaks droso
library(stringr)
library(data.table)
library(ROCR)
library(ggplot2)
library(dplyr)
library(patchwork)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/drosophila/modERN"
pathPoly="/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/polymorphism_analyses/NarrowPeaks/drosophila/modERN/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
sp = "drosophila"
minsub=5

method="exact_ranked_ancestral" 
AllMerged = paste0(path, "allMLE_drosophila_", method, "_", minsub, "sub_all_exp.Rds")
col=c("forestgreen", "orange", "deepskyblue4", "firebrick")
names(col) = c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model")

all_MLE_list <- readRDS(AllMerged)
all_MLE <- do.call(rbind, all_MLE_list)
all_MLE$peakID <- unlist(lapply(all_MLE_list, rownames))
rownames(all_MLE) <- all_MLE$peakID

# Filter Low Sample Size
all_MLE <- all_MLE[which(all_MLE$original_sample_size > 1500 & all_MLE$AUC >0.8),]

# Determine Conclusion according to FDR and alpha Tresholds
fdr_threshold=0.05
pval_threshold=0.01
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_threshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_threshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
all_MLE$Conclusion <- factor(all_MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model"))

all_MLE$signif <- as.factor(ifelse(all_MLE$pval.high <= pval_threshold, "Directional (+)", ifelse(all_MLE$pval.high >= 1-pval_threshold, "Directional (-)", "Neutral")))
all_MLE$signif <- factor(all_MLE$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$both <- grepl("Directional", all_MLE$signif) & grepl("Directional", all_MLE$FDR_Conclusion)

#### A - Proportion of peaks under selection
Conclu.sample <- prop.table(table(all_MLE$signif, all_MLE$exp), margin=2)
Conclu.sample.FDR <- prop.table(table(all_MLE$FDR_Conclusion, all_MLE$exp), margin=2)
Conclu.sample.both <- prop.table(table(all_MLE$both, all_MLE$exp), margin=2)

# Keep only Directional +
#Conclu.sample.FDR <- Conclu.sample.FDR[Conclu.sample.FDR$Var1 %in% c("Directional (+)"), ]

dfA <- data.frame(Proportion = c(Conclu.sample[1,] + Conclu.sample[2,], Conclu.sample.FDR[1,] + Conclu.sample.FDR[2,]),
  Method = rep(c("Permutations", "RegEvol"),  each = length(Conclu.sample[1,])))

pA <- ggplot(dfA, aes(x = Method, y = Proportion, fill = Method)) +
  geom_violin(alpha = 0.4, color = "black", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("RegEvol" = "red", "Permutations" = "black")) +
  labs(x = "", y = "Proportion Directional", fill = "") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )

pA
medDir = median(Conclu.sample.FDR[1,] + Conclu.sample.FDR[2,])
prop_Dir = nrow(all_MLE[grepl("Directional", all_MLE$FDR_Conclusion),])/nrow(all_MLE)
message("Median Proportion Directional RegEvol: ", medDir)
message("Proportion Directional RegEvol: ", prop_Dir)

#### B - Proportion of peaks under selection per ΔSVM quantile class
# Cut ΔSVM into quantile classes
breaks <- quantile(all_MLE$deltaSVM, probs = seq(0, 1, 0.05), na.rm = TRUE)
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2
all_MLE$classdeltaSVM <- cut(all_MLE$deltaSVM, breaks=breaks, labels = names(midpoints), include.lowest = T)

# Permutation-based proportions
perm_tab <- prop.table(table(all_MLE$signif, all_MLE$classdeltaSVM), margin = 2)
perm_df <- as.data.frame(perm_tab)
names(perm_df) <- c("Model", "DeltaClass", "Proportion")
perm_df$Method <- "Permutations"

# RegEvol-based proportions
reg_tab <- prop.table(table(all_MLE$Conclusion, all_MLE$classdeltaSVM), margin = 2)
reg_df <- as.data.frame(reg_tab)
names(reg_df) <- c("Model", "DeltaClass", "Proportion")
reg_df$Method <- "RegEvol"

# Combine
df_B <- rbind(perm_df, reg_df)
df_B <- df_B[df_B$Model == "Directional (+)", ]

# --- Plot B ---
x_labels <- levels(df_B$DeltaClass)
keep_labels <- rep("", length(x_labels))
keep_labels[c(1, seq(5, length(x_labels) - 1, by = 5), length(x_labels))] <- x_labels[c(1, seq(5, length(x_labels) - 1, by = 5), length(x_labels))]

pB <- ggplot(df_B, aes(x = DeltaClass, y = Proportion, color = Method, group=Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("RegEvol" = "red", "Permutations" = "black")) +
  labs(x = expression(Delta*"SVM (quantile)"), y = "Proportion Directional", color = "") +
  theme_minimal(base_size = 16) +
  scale_x_discrete(labels = keep_labels) +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 14),
    legend.position = c(0.2, 0.85),
    strip.text = element_text(size = 14, face = "bold"))


#### C & D - Polymorphism
TFS = levels(as.factor(all_MLE$exp))
fisher_df <- data.frame(TF = TFS, odds_ratio = NA, p_value = NA)

for (i in seq_along(TFS)){
  TF <- TFS[i]
  print(TF)
  sub_MLE <- all_MLE[which(all_MLE$exp==TF),]
  # Directional vs Other
  pos <- sub_MLE[grepl("Directional", sub_MLE$FDR_Conclusion),]
  nonPos <- sub_MLE[!grepl("Directional", sub_MLE$FDR_Conclusion),]
  subNumb <-c(sum(pos$Nmut),sum(nonPos$Nmut))
  polyNumb <-c(sum(pos$NbSNP),sum(nonPos$NbSNP))
  
  fisher_df$Npos[i] <- nrow(pos)
  fisher_df$NnoPos[i] <- nrow(nonPos)
  fisher_df$posMut[i] <- sum(pos$Nmut)
  fisher_df$posSNP[i] <- sum(pos$NbSNP)
  fisher_df$nonPosMut[i] <- sum(nonPos$Nmut)
  fisher_df$nonPosSNP[i] <- sum(nonPos$NbSNP)
  
  test <- fisher.test(matrix(c(subNumb,polyNumb),nrow = 2,ncol = 2))
  fisher_df$odds_ratio[i] <- test$estimate
  fisher_df$p_value[i] <- test$p.value
}
rownames(fisher_df) <- fisher_df$TF

# --- C: Barplot as ggplot ---
TF <- "CG8478_yellow_cinnabar_brown_speck_embryonic_1"
main <- "embryonic TwdlD peaks"

df_bar <- tibble(Category = c("Directional", "Other"),
  Sub_per_Poly = c(fisher_df[TF, "posMut"], fisher_df[TF, "nonPosMut"]) / c(fisher_df[TF, "posSNP"], fisher_df[TF, "nonPosSNP"]),
  n_label = c(fisher_df[TF, "Npos"], fisher_df[TF, "NnoPos"]))

pC <- ggplot(df_bar, aes(x = Category, y = Sub_per_Poly, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = paste0("n=", n_label)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Directional" = "red", "Other" = "deepskyblue3")) +
  ylim(0, 1.3) +  labs(y = "# Substitutions / # Polymorphisms", x = "") +
  ggtitle(main) + theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", margin = margin(b = 5)),  # small bottom margin
    plot.title.position = "panel",
    legend.position = "none",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) +
  annotate("text", x = 2.3, y = 1.25, size = 5, hjust = 1,
           label = paste0("ratio=", signif(fisher_df[TF, "odds_ratio"], 2),  "; p-val=", signif(fisher_df[TF, "p_value"], 2)))

# --- D: Volcano plot as ggplot ---
# Filter
fisher_df2 <- fisher_df %>% filter(odds_ratio > 0, Npos > 10)
fisher_df2$FDR <- p.adjust(fisher_df2$p_value, method = "BH")

# Create color factor
fisher_df2 <- fisher_df2 %>%
  mutate(color_group = case_when(
    FDR < 0.05 & odds_ratio > 1 ~ "Directional >",
    FDR < 0.05 & odds_ratio < 1 ~ "Non-directional >",
    TRUE ~ "NS"))

pD <- ggplot(fisher_df2, aes(x = odds_ratio, y = -log10(FDR), color = color_group)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Directional >" = "red",
                                "Non-directional >" = "deepskyblue3",
                                "NS" = "grey")) +
  labs(x = "Odds Ratio", y = expression(-log[10](FDR))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = c(0.25, 0.8))

# Display plots
combine <- (pA | pB) / (pC | pD) / (p1 | p2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
combine

# Save figure
ggsave(paste0(pathFigures, "Figure5_4sub_drosophila.png"), plot = combine, width = 11, height = 14, dpi = 320)
ggsave(paste0(pathFigures, "Figure5_4sub_drosophila.pdf"), plot = combine, width = 11, height = 14, dpi = 320)

