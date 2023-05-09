library(data.table)
library(reshape2)
library(ggplot2)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species <- c("human", "macaca", "mouse", "rat", "dog", "chicken", "cat", "rabbit")
TFs <- c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

models = read.table(paste0(path, "positive_selection/human/CEBPA_model_prediction.txt"), row.names = 1)
colnames(models) <- "human_CEBPA"
                    
for (sp in species){
  for (TF in TFs){
    file = paste0(path, "positive_selection/", sp, "/", TF, "_model_prediction.txt")
    if (file.exists(file)){
      model = read.table(file, row.names = 1)
      models[[paste0(sp, "_", TF)]] <- model$V2
    }
  }
}

cormat <- round(cor(models),2)

# Reorder the matrix: use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]

melted_cormat <- melt(cormat, na.rm = TRUE)

pdf(paste0(path, "figures/models_correlation.pdf"))
ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()

dev.off()

