library(ggplot2)
library(dplyr)
library(patchwork)

# Load Data
pathFigure <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
method="exact_ranked_ancestral" 

AllMerged = paste0(path, "drosophila_modERN_RegEvol_exact_ranked.Rds")
all_MLE_list <- readRDS(AllMerged)
all_MLE <- do.call(rbind, all_MLE_list)

original_data <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/modERN/FlyTFPeaksPrimaryTargets.tsv"))
gene_match <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/gene_match.txt"))
gene_names = fread(paste0(path, "../data/drosophila_gene_names.txt"), h=F)
colnames(gene_names) = c("Ensembl.Gene.ID", "Gene.name")

# Merge both Dataframes 
original_data$peakID <- original_data[["chrom:chromStart:chromEndID:experiment"]]
original_data <- original_data[,c("peakID", "Gene", "targetDist")]
fullData <- merge(all_MLE, original_data, by = "peakID")

# All genes
all_genes_original = as.data.frame(table(fullData$Gene))
colnames(all_genes_original) = c("Gene.name", "Freq")
rownames(all_genes_original) = all_genes_original$Gene.name
all_genes_original$MeanSub <- tapply(fullData$Nmut, fullData$Gene, mean)

# All TFs
all_TF_original = as.data.frame(table(fullData$TF))
colnames(all_TF_original) = c("Gene.name", "Freq")
rownames(all_TF_original) = all_TF_original$Gene.name

# Omega from Selectome
Selectome = fread("/Users/alaverre/Documents/Detecting_positive_selection/results/SQL_drosophila_omega.csv")
colnames(Selectome) = c("Ensembl.Gene.ID", "Gene.name_2", "omega0", "omega2", "lrt", "qvalue")
Selectome <- unique(merge(Selectome, gene_names, by = "Ensembl.Gene.ID"))
Selectome$lrt <- log(Selectome$lrt+1)

# Gene expression change
anat = "embryo"
gene.expr = fread(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/gene_expression/droso_", anat, "_stats.txt"))
gene.expr <- unique(merge(gene.expr, gene_names, by = "Ensembl.Gene.ID"))
gene.expr$abs_log2FC <- abs(gene.expr$log2FC)
fullData$meanlog2_mel <- gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$meanlog2_mel
fullData$abs_log2FC <- abs(gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$log2FC)
fullData$FDR <- -log10(gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$FDR)


evol_data <- Selectome # or gene.expr
stat <- "lrt" # or "abs_log2FC" or "meanlog2_mel"
genes <- "all" # or TF
################################################################################
# Merge and bin
if (genes == "all"){
  analysed_genes <- all_genes_original
  breaks = c(0, 10, 25, 100, 500, 1000, max(all_genes_original$Freq))
  labels = c("0-10", "11-25", "26-100", "101-500", "501-1000", ">1001")
}else{
  analysed_genes <- all_TF_original
  breaks = c(0, 2000, 3000, 5000, 7500, 10000, max(all_TF_original$Freq))
  labels = c("1-2k", "2-3k", "3-5k", "5-7.5k", "7.5-10k", ">10k")
}


all_genes <- merge(analysed_genes, evol_data, by = "Gene.name") %>%
  mutate(Freq_bin = cut(Freq, breaks = breaks, include.lowest = TRUE,
                        labels = labels))
# Spearman correlation
cor <- cor.test(all_genes$Freq, all_genes[[stat]], method = "spearman")

ylab = ifelse(stat == "meanlog2_mel", "Mean gene expression log2(TPM)", "log2 Fold Change") # 
ylab = stat

p1 <- ggplot(all_genes, aes(x = Freq_bin, y = !!sym(stat))) +
  geom_violin(fill = "gray70", color = NA, alpha = 0.6, trim = TRUE) +       
  geom_boxplot(outlier.shape = NA, notch = TRUE,
               fill = "white", alpha = 0.7, width = 0.5,                     
               color = "black", size = 0.6) +                               
  theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(x = "Nb peaks", y = ylab) +
  annotate("text", x=3.5, y = max(all_genes[[stat]], na.rm = TRUE) * 0.95,
           label = paste0("Spearman's rho = ", round(cor$estimate, 2),
                          "\n p-val = ", format.pval(cor$p.value, digits = 2)), size = 5, hjust = 0.5)


################################################################################
# Ratio Positive
# Get Positive according to Test
pos_genes = as.data.frame(table(fullData[fullData$FDR_Conclusion == "Directional (+)",]$Gene))
colnames(pos_genes) = c("Gene.name", "Freq")
rownames(pos_genes) = pos_genes$Gene.name
pos_genes$all <- all_genes[match(pos_genes$Gene.name, all_genes$Gene.name),]$Freq
pos_genes_filtered <- pos_genes %>% 
  filter(all >= 2) %>%
  mutate(ratio = Freq / all,
         ratio_bin = cut(ratio, breaks = c(0, 0.05, 0.1, 0.2, 1), include.lowest = TRUE))

# Spearman correlation
cor_res <- cor.test(pos_genes_filtered$ratio, pos_genes_filtered[[stat]], method = "spearman")

# Boxplot with points
p2 <- ggplot(pos_genes_filtered, aes(x = ratio_bin, y = !!sym(stat))) +
  geom_violin(fill = "gray70", color = NA, alpha = 0.6, trim = TRUE) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = "white", color = "black", alpha=0.7) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  labs(x = "Ratio Positive", y = "") +
  annotate("text", x = 2.5, y = max(pos_genes_filtered[[stat]], na.rm = TRUE) * 0.95,
    label = paste0("Spearman's rho = ", round(cor_res$estimate, 2), 
                   "\n p-val = ", format.pval(cor_res$p.value, digits=2)),
    size = 5,  hjust = 0.5)

# Combine plots
final_plot <- p1 + p2 
print(final_plot)
# Save plot
#ggsave(paste0(pathFigure, "SupFig_droso_", omega, "_", anat, ".png"), 
#       plot = final_plot, width = 14, height = 7, dpi = 300)

################################################################################