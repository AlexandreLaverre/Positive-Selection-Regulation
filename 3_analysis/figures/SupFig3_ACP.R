library(data.table)
library(progress)
library(factoextra)
library(patchwork)
library(Biostrings)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
species <- c("human", "macaca", "mouse", "spretus", "caroli", "rat", "dog", "chicken", "cat", "rabbit", "drosophila")
samples <- c("Wilson", "Schmidt12", "Ni12", "Rensch", "Myers", "Schmidt10", "Stefflova")
TFs <- c("HNF4A", "HNF6", "CEBPA", "CTCF", "FOXA1")

#species = c("drosophila")
#samples = c("modERN")
#TFs <- list.dirs(paste0(path, "drosophila/modERN"), recursive = FALSE)

get_gc_content <- function(kmers) {
  sapply(kmers, function(k) {
    gc <- sum(strsplit(k, "")[[1]] %in% c("G", "C"))
    return(gc / nchar(k))
  })
}

########################################################################
# 1 - Get all datas
AllModels <- paste0(path, "all_models.rds")

if (file.exists(AllModels)){
  # Load data
  models <- readRDS(AllModels)
} else {
  models = read.table(paste0(path, "human/Wilson/FOXA1/Model/kmer_predicted_weight.txt"), row.names = 1)
  colnames(models) <- "human   FOXA1"
  pb <- progress_bar$new(format = "  Processing [:bar] :percent (:elapsed s) | ETA :eta", total = length(TFs), clear = FALSE, width = 60)
  
  for (sp in species){
    for (TF in TFs){
      TF = basename(TF) # for drosophila
      for (sample in samples){
        file = paste(path, sp, sample, TF, "Model/kmer_predicted_weight.txt", sep="/")
        if (file.exists(file)){
          if (sample == "modERN"){
            pb$tick()
            simpleID <- unlist(strsplit(TF, "_"))
            simpleID <- paste(c(simpleID[1], tail(simpleID, 2)), collapse = "_")
          }else{
            message(paste("Opening:", sp, TF))
            simpleID <- paste(sp, TF, sep = " ")
          }

          model = read.table(file, row.names = 1)
          models[[simpleID]] <- model$V2
          #top <- DNAStringSet(head(rownames(model[order(model$V2, decreasing=T),, drop=F]),n=100))
          #writeXStringSet(top, paste(path, sp, sample, TF, "Model/top_kmer.fa", sep="/"), format = "fasta")
        }
      }
    }
  }
  models <- models[,-1]
  # Save data
  saveRDS(models, file=AllModels)
}

if ("modERN" %in% samples){
  exp_info <- read.csv("/Users/alaverre/Documents/Detecting_positive_selection/cluster/data/droso_exp_info.csv", h=T)
  base_names <- paste(exp_info$Gene, exp_info$Stage, sep = "_")
  name_counts <- ave(base_names, base_names, FUN = seq_along)
  unique_names <- paste(base_names, name_counts, sep = "_")
  rownames(exp_info) <- unique_names
  exp_info <- exp_info[colnames(models),]
}else{
  exp_info <- read.table("/Users/alaverre/Documents/Detecting_positive_selection/results/peaks_calling/peaks_numbers.txt", h=T)
  rownames(exp_info) <- paste(exp_info$species, exp_info$TF, sep = " ")
}

########################################################################
# 2 - PCA 
models_kmer <- t(models)
res.pca <- PCA(models_kmer, scale.unit = F)
#res.pca_filtered <- PCA(models_kmer_filtered, scale.unit = T)

names_PCA <- rownames(res.pca$ind$coord)
# Scree plot (variance explained)
fviz_eig(res.pca, addlabels = TRUE)

# PCA plot colored by TF (can be replaced by Antibody for drosophila)
make_pca_plot <- function(res.pca, axes, exp_info, names_PCA, panel_label) {
  fviz_pca_ind(
    res.pca, axes = axes,
    habillage = factor(exp_info[names_PCA, ]$TF),
    addEllipses = TRUE, ellipse.level = 0.95,
    palette = "npg",
    title = "") +
    theme_minimal(base_size = 16) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0)
    )
}

# Generate both plots
pA <- make_pca_plot(res.pca, c(1, 2), exp_info, names_PCA, "A") + guides(color = "none", fill = "none", shape = "none")
pB <- make_pca_plot(res.pca, c(2, 3), exp_info, names_PCA, "B")

########################################################################
# Correlate SVM with GC content overall
get_gc_content <- function(kmers) {
  sapply(kmers, function(k) {
    gc <- sum(strsplit(k, "")[[1]] %in% c("G", "C"))
    return(gc / nchar(k))
  })
}

mean_svm_per_kmer <- colMeans(models_kmer)
kmer_gc_ratio <- get_gc_content(rownames(models))

pc1_values <- res.pca$ind$coord[,1]
gc_weighted <- models_kmer %*% kmer_gc_ratio

cor_test <- cor.test(gc_weighted[names(pc1_values),], pc1_values, method="spearman")
rho <- round(cor_test$estimate, 3)
pval <- signif(cor_test$p.value, 3)

GC_data <- data.frame(PC1 = pc1_values, GC_weighted = gc_weighted[names(pc1_values),])
panelC <- ggplot(GC_data, aes(x = PC1, y = GC_weighted)) +
  geom_point() + theme_minimal(base_size = 14) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(x = "Dim1", y = "GC-weighted SVM score") +
  annotate("text", x = min(pc1_values), y = max(gc_weighted), 
           label = paste0("Spearman's rho = ", rho, "\np-values = ", pval), hjust = 0, vjust = 1)

########################################################################
# Save figure
combine <- ((pA | panelC) / pB) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18, face = "bold")) &
              plot_layout(guides = "collect") &  theme(legend.position = "bottom")

ggsave(paste0(pathFigures, "SupFigure3_PCA.png"), plot = combine, width = 10, height = 12, dpi = 320)
ggsave(paste0(pathFigures, "SupFigure3_PCA.png"), plot = combine, width = 10, height = 12, dpi = 320)

