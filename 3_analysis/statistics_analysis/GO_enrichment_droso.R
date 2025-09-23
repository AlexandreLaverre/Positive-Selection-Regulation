# Gene Ontoloy Enrichment
library(data.table)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)

# Load Data
path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
method="exact_ranked_ancestral" 
all_MLE <- fread(paste0(path, "positive_selection/NarrowPeaks/drosophila/allMLE_drosophila_", method, "_filtered.tsv"))
original_data <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/modERN/save/FlyTFPeaksPrimaryTargets.peaks.bed"))
gene_match <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/gene_match.txt"))
  
# Merge both Dataframes 
original_data$peakID <- paste(gsub("chr", "", original_data$chrom), original_data$chromStart, original_data$chromEnd, original_data$experiment, sep=":")
original_data <- original_data[,c("peakID", "Gene")]
fullData <- merge(all_MLE, original_data, by = "peakID")

all_genes_original = as.data.frame(table(fullData$Gene))
colnames(all_genes_original) = c("Gene.name", "Freq")
rownames(all_genes_original) = all_genes_original$Gene.name

# Mean mutations per gene
all_genes_original$MeanSub <- tapply(fullData$Prop.Sub, fullData$Gene, mean)

# Omega
omega = fread(paste0(path, "../data/drosophila_melanogaster_omega.txt"))
gene_names = fread(paste0(path, "../data/drosophila_gene_names.txt"), h=F)
colnames(gene_names) = c("Ensembl.Gene.ID", "Gene.name")
omega <- merge(omega, gene_names, by = "Ensembl.Gene.ID")
fullData$omega <- omega[match(fullData$Gene, omega$Gene.name),]$omega
fullData$lrt <- omega[match(fullData$Gene, omega$Gene.name),]$lrt

# New omega
Selectome = fread("/Users/alaverre/Documents/SQL_drosophila_omega.csv")
colnames(Selectome) = c("Ensembl.Gene.ID", "Gene.name", "omega0", "omega2", "lrt", "qvalue")
Selectome <- unique(merge(Selectome, gene_names, by = "Ensembl.Gene.ID"))
Selectome$Gene.name <- Selectome$Gene.name.y
Selectome$Gene.name.x <- NULL
Selectome$Gene.name.y <- NULL
fullData$omega0 <- Selectome[match(fullData$Gene, Selectome$Gene.name),]$omega0
fullData$omega2 <- Selectome[match(fullData$Gene, Selectome$Gene.name),]$omega2
fullData$lrt <- Selectome[match(fullData$Gene, Selectome$Gene.name),]$lrt


evol_data <- Selectome
omega <- "omega2"
pdf(paste0(path, "figures/drosophila_evol_", omega, "_Selectome.pdf"), width = 10, height = 6)
# Correl nb total and omega
all_genes <- merge(all_genes_original, evol_data, by = "Gene.name")
boxplot(all_genes[[omega]] ~ cut(all_genes$Freq, breaks = c(0, 10, 25, 100, 500, 1000)), notch=T, outline=F,
        xlab = "Nb Peaks", ylab = "Omega", main = "Nb Total Peaks vs Omega", cex=0.3)

cor <- cor.test(all_genes$Freq, all_genes[[omega]], method = "spearman")
mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)

# Correl nb substitution and omega
boxplot(all_genes[[omega]] ~ cut(all_genes$MeanSub, breaks = c(0, 0.025, 0.05, 0.075, 0.1, 0.2)), notch=T, outline=F,
        xlab = "Nb Peaks", ylab = omega, main = "Nb Total Peaks vs Omega", cex=0.3)

cor <- cor.test(all_genes$Freq, all_genes[[omega]], method = "spearman")
mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)

# Ranked GO - Nb Total
geneList <- all_genes$Freq
names(geneList) <- all_genes$Gene.name
geneList <- sort(geneList, decreasing = TRUE)

#gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
#               pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 1, 
#               maxGSSize = 500, scoreType = "pos", eps = 0)
#dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Nb Total")

#test <- tapply(fullData$FDR_null_pos, fullData$Gene, function(x) which(p.adjust(x, method = "bonferroni")<0.5))
######################## Get Gene Lists ######################## 
genes_ID <- list()
genes_ID_background <- list()
#par(mfrow=c(2,1))
for (test in c("Conclusion", "FDR_Conclusion", "signif", "signif.FDR")){
  print(test)
  # Get Positive according to Test
  pos_genes = as.data.frame(table(fullData[get(test) == "Directional (+)",]$Gene))
  colnames(pos_genes) = c("Gene.name", "Freq")
  rownames(pos_genes) = pos_genes$Gene.name

  pos_genes <- merge(pos_genes, evol_data, by = "Gene.name")
  boxplot(pos_genes[[omega]] ~ cut(pos_genes$Freq, breaks = c(0, 2, 10, 25, 100, 500)), notch=T, outline=F,
          xlab = "Nb Positive Peaks", ylab = omega, main=test, cex=0.3)
  
  cor <- cor.test(pos_genes$Freq, pos_genes[[omega]], method = "spearman")
  mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)
  
  # Nb Positive - Classic GO
  go <- enrichGO(gene = as.character(pos_genes$Gene.name), universe =  as.character(all_genes$Gene.name),
                 OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  
  #p1 <- dotplot(go, showCategory = 20) + ggtitle("Classic GO (>= 1 Positive)")
  
  # Ranked GO 
  geneList <- pos_genes$Freq
  names(geneList) <- pos_genes$Gene.name
  geneList <- sort(geneList, decreasing = TRUE)
  
  #gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
  #               pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 1, 
  #               maxGSSize = 500, scoreType = "pos", eps = 0)
  #p2 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Nb Positive")
  
  #emapplot(pairwise_termsim(gsego))
  # Ranked GO
  #go_r <- GSEA(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
  #             pvalueCutoff = 0.05, minGSSize = 1, maxGSSize = 500, pAdjustMethod = "BH")
  
  
  # Ranked GO - Ratio of positive / total
  pos_genes$all = all_genes[match(pos_genes$Gene.name, all_genes$Gene.name),]$Freq
  pos_genes_filtered <- pos_genes[pos_genes$all >= 2,]
  pos_genes_filtered$ratio = pos_genes_filtered$Freq / pos_genes_filtered$all
  boxplot(pos_genes_filtered[[omega]] ~ cut(pos_genes_filtered$ratio, breaks = c(0, 0.05, 0.1, 0.2, 1)), notch=T, 
          xlab = "Ratio Positive", ylab = omega, main=test, cex=0.3)
  
  cor <- cor.test(pos_genes_filtered$ratio, pos_genes_filtered[[omega]], method = "spearman")
  mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)
  
  geneList <- pos_genes_filtered$ratio
  names(geneList) <- pos_genes_filtered$Gene.name
  geneList <- sort(geneList, decreasing = TRUE)
  
  #gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
  #               pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = 1, 
  #               maxGSSize = 500, scoreType = "pos", eps = 0)
  #p3 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ratio positive/total")
  
  genes <- pos_genes_filtered[pos_genes_filtered$ratio >= 0.25,]$Gene.name
  genes_ID[[test]] <- gene_match[gene_match$Gene.name %in% genes,]$Gene.Stable.ID
  genes_ID_background[[test]] <- gene_match[gene_match$Gene.name %in% pos_genes_filtered$Gene.name,]$Gene.Stable.ID
  # write genes_ID
  #write.table(genes_ID, paste0(path, "figures/drosophila_genes_ratio_0.25_", test, ".txt"), 
  #            sep = "\t", row.names = F, col.names = F, quote = F)
  #write.table(genes_ID_background, paste0(path, "figures/drosophila_genes_ratio_", test, ".txt"), 
  #           sep = "\t", row.names = F, col.names = F, quote = F)
  
  # Plot all figures in the same pannel
  #pdf(paste0(path, "figures/drosophila_GO_enrichment_", test, ".pdf"), width = 14, height = 14)
  #plot <- (p1 | p2) / (p3)
  #print(plot) 
  #dev.off()
  
}
dev.off()

######################## Top Anat ######################## 
library(BgeeDB)

Enrich_classic <- list()
Enrich_elim <- list()
Enrich_weight <- list()

# Create Bgee object for Drosophila
bgee <- Bgee$new(species = "7227", dataType = "rna_seq")
myTopAnatData <- loadTopAnatData(bgee)

for (test in c("Conclusion", "FDR_Conclusion", "signif", "signif.FDR")){
  print(test)
  # Create gene list (1 = of interest, 0 = background)
  geneList <- factor(as.integer(c(genes_ID[[test]], genes_ID_background[[test]]) %in% genes_ID[[test]]))
  names(geneList) <- c(genes_ID[[test]], genes_ID_background[[test]])
  
  # Build TopAnat object
  myTopAnatObject <- topAnat(myTopAnatData, geneList)
  
  # Run different algorithms
  results_classic <- runTest(myTopAnatObject, algorithm = 'classic', statistic = 'fisher')
  results_elim    <- runTest(myTopAnatObject, algorithm = 'elim', statistic = 'fisher')
  results_weight  <- runTest(myTopAnatObject, algorithm = 'weight01', statistic = 'fisher')
  
  # Create results tables
  table_classic <- makeTable(myTopAnatData, myTopAnatObject, results_classic, cutoff = 1)
  table_elim    <- makeTable(myTopAnatData, myTopAnatObject, results_elim, cutoff = 1)
  table_weight  <- makeTable(myTopAnatData, myTopAnatObject, results_weight, cutoff = 0.5)

  # save
  if (length(table_classic) > 1){
    Enrich_classic[[test]] <- table_classic
  }
  if (length(table_elim) >1){
    Enrich_elim[[test]] <- table_elim
    Enrich_weight[[test]] <- table_weight # [,c(2,6,8)]
  }
  
}

######################## Plot TopAnat results ######################## 
library(ggplot2)
library(dplyr)

# Keep significant results only
top_terms <- Enrich_weight[[test]] %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(20)  # top 20 enriched terms

# Barplot
ggplot(top_terms, aes(x = reorder(organName, -log10(FDR)), y = -log10(FDR))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Anatomical Term", y = "-log10(FDR)", title = "TopAnat Enrichment (elim + Fisher)") +
  theme_minimal()

ggplot(top_terms, aes(x = foldEnrichment, y = reorder(organName, -log10(pValue)))) +
  geom_point(aes(size = annotated, color = -log10(FDR))) +
  scale_color_gradient(low = "lightblue", high = "red") +
  labs(x = "Fold Enrichment", y = "Anatomical Term",
       title = "TopAnat Dotplot (elim)",
       color = "-log10(FDR)", size = "Total Genes in Term") +
  theme_minimal()

