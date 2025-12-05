library(data.table)
library(reshape2)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(Biostrings)
library(progress)
library(pheatmap)
library(umap)
library(dplyr)
library(patchwork)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/NarrowPeaks/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
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
# 2 - Data processing
models_kmer <- t(models)

# Remove lowly variable kmer
var <- apply(models_kmer, 2, sd)
models_kmer_filtered <- models_kmer[, var > 0.25]

# Z-score scaling 
models_kmer_scaled <- scale(models_kmer_filtered)

# Transpose back
models_sample_filtered <- t(models_kmer_filtered)
models_sample_scaled <- t(models_kmer_scaled)


########################################################################
# Correlate SVM with GC content overall
mean_svm_per_kmer <- colMeans(models_kmer)
kmer_gc_ratio <- get_gc_content(rownames(models))

cor(mean_svm_per_kmer, kmer_gc_ratio, method="spearman")
kmer_gc_ratio_bin = cut(kmer_gc_ratio, breaks = seq(0, 1, 0.05), include.lowest = TRUE)
boxplot(mean_svm_per_kmer~kmer_gc_ratio_bin, outline=F, notch=T, xlab="GC content", ylab="Mean SVM score")

# By TF
TF_names <- sub(".* ", "", rownames(models_kmer))
plot_list <- list()
cor_list <- list()

for (tf in unique(TFs)) {
  rows <- which(TF_names == tf)
  mean_svm_per_kmer <- colMeans(models_kmer[rows, , drop=FALSE])
  
  df_tf <- data.frame(TF = tf, GC_bin = kmer_gc_ratio_bin, GC = kmer_gc_ratio, SVM = mean_svm_per_kmer)
  
  # Spearman correlation
  cor_val <- cor(mean_svm_per_kmer, kmer_gc_ratio, method = "spearman")
  cor_list[[tf]] <- data.frame(TF=tf, cor=round(cor_val, 3))
  
  plot_list[[tf]] <- df_tf
}

plot_df <- bind_rows(plot_list)
cor_df <- bind_rows(cor_list)

# Boxplot with facet per TF
p <- ggplot(plot_df, aes(x = GC_bin, y = SVM)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  facet_wrap(~TF, scales = "fixed") +   # same scale across facets
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "GC content bin", y = "Mean SVM score") +
  coord_cartesian(ylim = c(-1.5, 0.5))      # fix y-axis to [-1, 1]

# Add correlation text
p + geom_text(
  data = cor_df,
  aes(x = 10, y = 0.5, label = paste0("rho=", cor)),
  inherit.aes = FALSE,
  hjust = 1
)

# Correlate PC1 with GC content
pc1_values <- res.pca$ind$coord[,1]
gc_weighted <- models_kmer %*% kmer_gc_ratio

cor_test <- cor.test(gc_weighted[names(pc1_values),], pc1_values)
rho <- round(cor_test$estimate, 3)
pval <- signif(cor_test$p.value, 3)

plot(pc1_values, gc_weighted[names(pc1_values),], pch=19, xlab="PC1", ylab="Weighted GC score")
abline(lm(gc_weighted[names(pc1_values),] ~ pc1_values), col="red")
legend("topleft", legend = paste0("rho = ", rho, "\nP = ", pval), bty = "n", text.col = "black")

pca_df <- data.frame(PC1 = res.pca$ind$coord[,1], PC2 = res.pca$ind$coord[,2], GC_weigth = gc_weighted[names_PCA,])
ggplot(pca_df, aes(x=PC1, y=PC2, color=GC_weigth)) +
  geom_point(size=2) +
  scale_color_gradient(low="blue", high="red") +
  theme_minimal()

########################################################################
# 5 - UMAP for droso
umap_result <- umap(models_kmer_scaled, n_components=2)

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = exp_info[rownames(res.pca$ind$coord),]$species)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP of k-mer SVM score", color = "z")

### Define clusters
clusters <- kmeans(umap_result$layout, centers = 3)$cluster
umap_df$Cluster <- factor(clusters)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2) + theme_minimal()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = exp_info$Peaks)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  theme_minimal() +
  labs(color = "Peak count",
       title = "UMAP colored by number of Peaks")


# Compare average score of k-mers per cluster
#cluster_kmer_avg <- aggregate(models_scaled, by = list(cluster = clusters), FUN = mean)

# Correlate k-mers with UMAP2
correlations <- apply(models_kmer_scaled, 2, function(x) cor(x, umap_result$layout[,2]))
cor_df <- data.frame(kmer = names(correlations), cor = correlations)
cor_df <- cor_df[order(abs(cor_df$cor), decreasing = TRUE), ]
head(cor_df, 20)  # Top k-mers correlated with UMAP2

kmer_of_interest <- "GCGCGTACTA"
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = -models_kmer_scaled[, kmer_of_interest])) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(color = kmer_of_interest)


########################################################################
# Plot all PCAs
pdf(paste0(pathFigures, "models_PCAs_all_sp.pdf"), onefile=T)
p <- fviz_pca_ind(res.pca, habillage = as.factor(rownames(models_t)), 
                  addEllipses=T, ellipse.level=0.95) + scale_color_brewer(palette="Dark2") + theme_minimal() + theme(legend.position = "none")
print(p)


dev.off()
########################################################################

fviz_eig(res.pca.scale, addlabels = TRUE, ylim = c(0, 80))  # from factoextra
fviz_pca_ind(res.pca,  geom = "point", repel = TRUE,
             col.ind = "cos2", # color by quality of representation
             gradient.cols = c("white", "blue", "red"))

fviz_pca_ind(res.pca, geom = "point", addEllipses = TRUE,
             habillage = as.factor(exp_info$Stage)) # factor of group/condition

pca_coords <- as.data.frame(res.pca$ind$coord)

ggplot(pca_coords, aes(x = Dim.2, y = Dim.3, color = as.factor(exp_info$Stage))) +
  geom_point(size = 2) +
  labs(x = "PC1", y = "PC2", color = "Group") +
  theme_minimal()

fviz_pca_var(res.pca, select.var = list(contrib = 1000))
fviz_pca_var(res.pca, col.var = "contrib") +
  scale_color_gradient2(low = "white", high = "red", midpoint = median(res.pca$var$contrib[,1]))


var <- apply(models_t, 2, sd)
mean <- apply(models_t, 2, mean)



HNF4A=datas[["all"]][["human    HNF4A"]]
HNF6=datas[["all"]][["human    HNF6"]]
CEBPA=datas[["all"]][["human    CEBPA"]]
CTCF=datas[["all"]][["human    CTCF"]]
FOXA1 = datas[["all"]][["human   FOXA1"]]

par(mfrow=c(2,2))
plot(FOXA1, CTCF, cex=0.1)
segments(x0=-10, y0=-10, x1=10, y1=10, col="red")
plot(FOXA1, HNF6, cex=0.1)
segments(x0=-10, y0=-10, x1=10, y1=10, col="red")
plot(FOXA1, CEBPA, cex=0.1)
segments(x0=-10, y0=-10, x1=10, y1=10, col="red")
plot(FOXA1, HNF4A, cex=0.1)
segments(x0=-10, y0=-10, x1=10, y1=10, col="red")

data = datas[["all"]]
top_HNF4A=rownames(data[order(data$`human    HNF4A`, decreasing = T),])[1:100]
top_FOXA1=rownames(data[order(data$`human    FOXA1`, decreasing = T),])[1:100]
top_common=top_FOXA1[which(top_FOXA1%in%top_HNF4A)]
names(top_common) <- top_common
top_common <- DNAStringSet(top_common, use.names=T)
writeXStringSet(top_common, paste(path, "/human/Wilson/HNF4A/Model/top_common_HNF4A_FOXA1.fa", sep="/"), format = "fasta")

high_common=rownames(data[which(data$`human   FOXA1`>2.5 & data$`human    HNF4A`>2.5),])
names(high_common) <- high_common
high_common=DNAStringSet(high_common, use.names=T)
writeXStringSet(high_common, paste(path, "/human/Wilson/HNF4A/Model/high_common_HNF4A_FOXA1.fa", sep="/"), format = "fasta")
