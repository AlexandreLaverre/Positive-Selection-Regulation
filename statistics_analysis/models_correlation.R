library(data.table)
library(reshape2)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(Biostrings)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species <- c("human", "macaca", "mouse", "spretus", "caroli", "rat", "dog", "chicken", "cat", "rabbit")

########################################################################
# Get all datas
datas.names <- c("all", "noCTCF", "noHNF6", "noCEBPA")
TFs <- c("FOXA1", "HNF4A", "CEBPA", "HNF6", "CTCF")
TFs_noCTCF <- c("FOXA1", "HNF4A", "CEBPA", "HNF6")
TFs_noHNF6 <- c("FOXA1", "HNF4A", "CEBPA")
TFs_noCEBPA <- c("FOXA1", "HNF4A")
TFs_list <- list(TFs, TFs_noCTCF, TFs_noHNF6, TFs_noCEBPA)
names(TFs_list) <- datas.names

datas <- vector("list", length(datas.names))
names(datas) <- datas.names
TF_class <- vector("list", length(datas.names))
names(TF_class) <- datas.names

for (TF_data in names(TFs_list)){
  models = read.table(paste0(path, "positive_selection/human/FOXA1_model_prediction.txt"), row.names = 1)
  colnames(models) <- "human   FOXA1"
  TF_var <- c()    
  for (sp in species){
    for (TF in TFs_list[[TF_data]]){
      file = paste0(path, "positive_selection/", sp, "/", TF, "_model_prediction.txt")
      if (file.exists(file)){
        model = read.table(file, row.names = 1)
        models[[paste0(sp, "    ", TF)]] <- model$V2
        TF_var <- c(TF_var, TF)
        
        top <- DNAStringSet(head(rownames(model[order(model$V2, decreasing=T),, drop=F]),n=25))
        writeXStringSet(top, paste0(path, "positive_selection/", sp, "/", TF, "_top_kmer.fa"), format = "fasta")
      }
    }
  }
  datas[[TF_data]] <- models
  TF_class[[TF_data]] <- TF_var
}


########################################################################
# Run all PCAs
PCAs <- vector("list", length(datas.names))
names(PCAs) <- datas.names

for (data in names(datas)){
  models_t <- t(datas[[data]])
  res.pca <- PCA(models_t, graph = T)
  PCAs[[data]] <- res.pca
}
#eig.val <- get_eigenvalue(res.pca)
#eig.val
#fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

#var <- get_pca_var(res.pca)
#var

########################################################################
# Models correlations clustering
cormat <- round(cor(datas[["all"]]), 2)

# Reorder the matrix: use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]

melted_cormat <- melt(cormat, na.rm = TRUE)

# Plot Heatmap
pdf(paste0(path, "figures/models_heatmap_correlation_mice_added.pdf"), width=9, height=9)
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1))+
  coord_fixed()
dev.off()

########################################################################
# Plot all PCAs
pdf(paste0(path, "figures/models_PCAs_mice_added.pdf"), onefile=T)
for (data in names(datas)){
  p <- fviz_pca_ind(PCAs[[data]], habillage = as.factor(TF_class[[data]]), 
                    addEllipses=T, ellipse.level=0.95) + scale_color_brewer(palette="Dark2") + theme_minimal()
  print(p)
}

dev.off()
########################################################################

