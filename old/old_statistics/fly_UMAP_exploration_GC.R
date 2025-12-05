# 5 - UMAP 
umap_result <- umap(models_kmer_scaled, n_components=2)

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = exp_info[rownames(res.pca$ind$coord),]$TF)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP of k-mer SVM score", color = "z")

### Define clusters
clusters <- kmeans(umap_result$layout, centers = 3)$cluster
umap_df$Cluster <- factor(clusters)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2) + theme_minimal()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = exp_info[rownames(res.pca$ind$coord),]$total_peaks)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  theme_minimal() +
  labs(color = "Peak count",
       title = "UMAP colored by number of Peaks")


# Simple function to get GC proportion in a k-mer
get_gc_content <- function(kmers) {
  sapply(kmers, function(k) {
    gc <- sum(strsplit(k, "")[[1]] %in% c("G", "C"))
    return(gc / nchar(k))
  })
}

# Apply it
gc_content <- get_gc_content(colnames(models_kmer_scaled))
mean_score <- apply(models_kmer_scaled, 2, mean)

boxplot(mean_score~gc_content, notch=T, outline=F,
        xlab="GC content of k-mer", ylab="Mean association score")


# Weighted GC content per dataset
gc_weighted <- models_kmer_scaled %*% gc_content
umap_df$gc_weighted <- as.vector(gc_weighted) 

# Visualize on UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = gc_weighted)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(color = "GC weight", title = "GC weigth gradient on UMAP")


# Split by UMAP2 median
top_samples <- umap_df$UMAP2 > median(umap_df$UMAP2)
bottom_samples <- !top_samples

avg_top <- colMeans(models_kmer_scaled[top_samples, ])
avg_bottom <- colMeans(models_kmer_scaled[bottom_samples, ])

kmer_diff <- avg_top - avg_bottom # Diff score
log2fc <- log2(avg_top + 1e-6) - log2(avg_bottom + 1e-6) # Log2 fold-change

top_kmers <- sort(kmer_diff, decreasing = T)[1:20]
bottom_kmers <- sort(kmer_diff, decreasing = F)[1:20]

par(mar=c(8,4,3,1))
barplot(top_kmers, las=2, col="#377eb8", main="Top 20 - Top UMAP2", ylab="Mean(Top) - Mean(Bottom)")
barplot(bottom_kmers, las=2, col="#e41a1c", main="Top 20 - Bottom UMAP2", ylab="Mean(Top) - Mean(Bottom)")

gc_of_top <- get_gc_content(names(top_kmers))
gc_of_bottom <- get_gc_content(names(bottom_kmers))

boxplot(gc_of_top, gc_of_bottom,
        names = c("Top UMAP2 (GC)", "Bottom UMAP2 (AT)"),
        main = "GC content of the most divergent kmers",
        ylab = "GC content", col = c("#377eb8", "#e41a1c"))

# Plot the top kmers
fc_df <- data.frame(kmer = names(kmer_diff), gc_content=gc_content, avg_top=avg_top, avg_bottom=avg_bottom, diff=kmer_diff, log2FC = log2fc)
fc_df <- fc_df[order(abs(fc_df$diff), decreasing = TRUE), ]
head(fc_df, 20)


top_kmers <- head(fc_df$kmer, 50)
mat <- models_t[, top_kmers]
pheatmap(mat, scale = "row", show_rownames = FALSE)
