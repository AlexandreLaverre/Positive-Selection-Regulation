library(stringr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Biostrings)
library(qvalue)
library(gap)
library(patchwork)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"

species=c("caroli", "spretus", "dog", "cat", "human", "macaca", "rabbit", "rat", "mouse", "chicken", "drosophila")
sample=c("Wilson", "Myers", "Rensch", "Schmidt10", "Schmidt12", "Stefflova", "Ni12")
TFs=c("CEBPA", "FOXA1", "HNF4A", "HNF6", "CTCF")
maxSub=150
minSub=4
pval_threshold=0.01
fdr_threshold=0.1

method="exact"
################################################################################
# Get an object of all merged results

AllMerged = paste0(path, "allMLE_list_", method, "_", minSub, "sub.Rds")
if (!file.exists(AllMerged)){
  all_MLE_list <- list()
  for (sp in species){
    for (samp in sample){
      for (TF in TFs){
        MLE.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/MLE_summary_", method, ".csv")
        delta.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/deltas/ancestral_to_observed_deltaSVM.txt")
        permut.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
        stats.file <- paste0(path, "positive_selection/NarrowPeaks/", sp, "/", samp, "/", TF, "/sequences/focal_substitutions_stats.txt")
        
        if (!file.exists(MLE.file)){next}
        print(paste(sp, samp, TF))
        
        MLE <- read.csv(MLE.file, h=T, row.names = 1)
        
        # Stats substitutions
        seq_stats <- read.table(stats.file, h=T, row.names = 1)
        seq_stats <- seq_stats[row.names(MLE),]
        
        # Observed deltas
        obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
        deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
        row.names(deltas) <- deltas$seq_name
        deltas <- deltas[row.names(MLE),]
        
        MLE$SVM <- deltas$SVM
        MLE$deltaSVM <- deltas$deltaSVM
        names <- gsub("_.*:\\d+:\\d+:", "_", rownames(MLE))
        MLE$ID <- gsub("[_\\:]([[:alpha:]]).*_.*", "", names)
        MLE$peaks_ID <- gsub(".*:\\d+:\\d+[_\\:]", "", rownames(MLE))
        
        # Filtered out according to indel and Nsub
        MLE$original_length <- as.numeric(str_split_i(MLE$ID, ":", 3))-as.numeric(str_split_i(MLE$ID, ":", 2))
        MLE$indel <- MLE$original_length - seq_stats$Length
        
        MLE$Diffdelta <- MLE$deltaSVM-MLE$SumObs
        MLE$Prop.indel <- (MLE$indel*100)/MLE$original_length
        MLE$original_sample_size <- nrow(MLE)
        MLE <- MLE[which(MLE$Prop.indel<10),]
        MLE <- MLE[which(MLE$Nmut >= minSub),]
        
        # FDR corrections
        MLE$FDR_null_pos <- p.adjust(MLE$p_value_null_pos, method="fdr")
        MLE$FDR_null_purif <- p.adjust(MLE$p_value_null_purif, method="fdr")
        MLE$FDR_purif_pos <- p.adjust(MLE$p_value_purif_pos, method="fdr")
        
        # Permut Test 
        permut <- read.table(permut.file, h=T, row.names = 1)
        permut <- permut[match(rownames(MLE), rownames(permut)),]
        MLE$permut.pval <- permut$pval.high
        MLE$permut.FDR <- p.adjust(permut$pval.high, method="fdr")
        
        # Annotations
        MLE$species <- sp
        if (samp=="Schmidt10"){samp="Schmidt12"}
        MLE$TF <- TF
        MLE$sample <- samp
        MLE$sampleSize <- nrow(MLE)
        
        rownames(MLE) <- MLE$ID
        all_MLE_list[[paste0(sp, "_", TF)]] <- MLE
      }
    }
  }
  saveRDS(all_MLE_list, file=AllMerged)
}else{all_MLE_list <- readRDS(AllMerged)
}
################################################################################
all_MLE <- do.call(rbind, all_MLE_list)

all_MLE$length <- as.numeric(str_split_i(all_MLE$ID, ":", 3))-as.numeric(str_split_i(all_MLE$ID, ":", 2))
all_MLE$Prop.Sub <- all_MLE$Nmut/all_MLE$length
species_order <- c("drosophila", "chicken", "cat", "rabbit", "rat", 
                   "spretus", "caroli", "mouse", "human", "macaca", "dog")

all_MLE$species <- factor(all_MLE$species, levels = species_order)

all_MLE$exp <- paste0(all_MLE$species, "_", all_MLE$TF)
all_MLE <- all_MLE %>% mutate(species = factor(species), TF = factor(TF) )
all_MLE$Group <- interaction(all_MLE$species, all_MLE$TF, sep = " - ")


################################################################################
# Determine Conclusion according to FDR and alpha Tresholds
all_MLE$FDR_Conclusion <- "Neutral"
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_null_purif <= fdr_threshold, "Stabilizing", all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- ifelse(all_MLE$FDR_purif_pos <= fdr_threshold, ifelse(all_MLE$AlphaPos>all_MLE$BetaPos,"Directional (+)", "Directional (-)"), all_MLE$FDR_Conclusion)
all_MLE$FDR_Conclusion <- factor(all_MLE$FDR_Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))
all_MLE$Conclusion <- factor(all_MLE$Conclusion, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral model", "Disruptive"))

all_MLE$signif <- as.factor(ifelse(all_MLE$permut.pval <= pval_threshold, "Directional (+)", ifelse(all_MLE$permut.pval >= 1-pval_threshold, "Directional (-)", "Neutral")))
all_MLE$signif <- factor(all_MLE$signif, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

all_MLE$signif.FDR <- "Neutral"
all_MLE$signif.FDR <- as.factor(ifelse(all_MLE$permut.FDR > fdr_threshold, "Neutral", ifelse(all_MLE$deltaSVM > 0, "Directional (+)", "Directional (-)")))
all_MLE$signif.FDR <- factor(all_MLE$signif.FDR, levels=c("Directional (+)", "Directional (-)", "Stabilizing", "Neutral"))

################################################################################
tf_palette <- ggpubr::get_palette("npg", 5)

# ---- A) Substitutions ----
pA <- ggplot(all_MLE, aes(x = Nmut, y = species, fill = TF)) +
  geom_boxplot(width = 0.9, color = "black", outlier.shape = NA, notch=T,
               position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = tf_palette, name = "TF") +
  labs(x = "Substitutions", y = "Species") +
  coord_cartesian(xlim = c(0, 50)) +   # <<< force max x to 50
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.border = element_blank()
  )

# ---- B) original_length----
pB <- ggplot(all_MLE, aes(x = original_length, y = species, fill = TF)) +
  geom_boxplot(width = 0.9, color = "black",  outlier.shape = NA, notch=T, 
               position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = tf_palette, name = "TF") +
  labs(x = "Peak length (bp)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right",
    panel.border = element_blank()
  )

# ---- C) SVM  ----
pC <- ggplot(all_MLE, aes(x = SVM, y = species, fill = TF)) +
  geom_boxplot(width = 0.9, color = "black", outlier.shape = NA, notch = TRUE,
               position = position_dodge2(preserve = "single") ) +
  scale_fill_manual(values = tf_palette, name = "TF") +
  labs(x = "SVM score", y = "Species") +
  xlim(-700, 300) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none", panel.border = element_blank())

# ---- D) deltaSVM  ----
pD <- ggplot(all_MLE, aes(x = deltaSVM, y = species, fill = TF)) +
  geom_boxplot(width = 0.9, color = "black", outlier.shape = NA, notch=T,
               position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = tf_palette, name = "TF") +
  labs(x = expression(Delta*"SVM"), y = NULL) +
  xlim(-25, 35) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.border = element_blank()
  )

# ---- Combine panels ----
combined <- ((pA | pB) / (pC | pD)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

combined

# Save figure
ggsave(paste0(pathFigures, "SupFig11_", minSub, "_sub_mammals.png"), plot = combined, width = 12, height = 14, dpi = 320)
ggsave(paste0(pathFigures, "SupFig11_", minSub, "_sub_mammals.pdf"), plot = combined, width = 12, height = 14, dpi = 320)

############################
# remove drosophila for proportion plot
all_MLE <- all_MLE %>% filter(species != "drosophila")

df_dir <- all_MLE %>% group_by(species, TF) %>%
  summarise(RegEvol = sum(signif.FDR %in% c("Directional (+)", "Directional (-)")),
            total_peaks = n(),
            prop_dir = RegEvol / total_peaks ) %>% ungroup()

# A - Proportion directional peaks per species and TF
pA <- ggplot(df_dir, aes(x = species, y = prop_dir, fill = TF)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = RegEvol), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  scale_fill_manual(values = tf_palette, name = "TF") +
  scale_y_continuous(limits = c(0, max(df_dir$prop_dir) * 1.1)) +
  labs(y = "Proportion of Directional peaks", x = "Species", fill = "TF") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=14),
    axis.text.y = element_text(size=14),
    legend.position = "bottom"
  )

# B - Proportion versus substitutions
df_dir_nmut <- all_MLE %>%
  mutate(Nmut_group = ifelse(Nmut > 30, "30+", as.character(Nmut))) %>%
  mutate(Nmut_group = factor(Nmut_group, 
                             levels = c(as.character(sort(unique(Nmut[Nmut <= 30]))), "30+"))) %>%
  group_by(Nmut_group, TF) %>%
  summarise(RegEvol = sum(signif.FDR %in% c("Directional (+)", "Directional (-)")),
            total_peaks = n(),
            prop_dir = RegEvol / total_peaks,
            .groups = "drop")

pB <- ggplot(df_dir_nmut, aes(x = Nmut_group, y = prop_dir, color = TF, group = TF)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = tf_palette, name = "TF") +
  labs(y = "Proportion of Directional peaks",
       x = "Number of substitutions") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14),
        axis.text.y = element_text(size=14),
        legend.position = "none")

combined <- (pA / pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

combined

ggsave(paste0(pathFigures, "SupFig12_", minSub, "_sub_mammals.png"), plot = combined, width = 10, height = 8, dpi = 320)
ggsave(paste0(pathFigures, "SupFig12_", minSub, "_sub_mammals.pdf"), plot = combined, width = 10, height = 8, dpi = 320)
