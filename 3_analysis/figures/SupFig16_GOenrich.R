# Gene Ontoloy Enrichment
library(data.table)
library(clusterProfiler)
library(org.Dm.eg.db)
library(enrichplot)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load Data
path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/"
method="exact_ranked_ancestral" 
all_MLE <- fread(paste0(path, "positive_selection/NarrowPeaks/drosophila/allMLE_drosophila_", method, "_2sub_all_exp.tsv"))
original_data <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/modERN/FlyTFPeaksPrimaryTargets.tsv"))
gene_match <- fread(paste0(path, "/peaks_calling/NarrowPeaks/drosophila/gene_match.txt"))
gene_names = fread(paste0(path, "../data/drosophila_gene_names.txt"), h=F)
colnames(gene_names) = c("Ensembl.Gene.ID", "Gene.name")

# Merge both Dataframes 
original_data$peakID <- original_data[["chrom:chromStart:chromEndID:experiment"]]
original_data <- original_data[,c("peakID", "Gene", "targetDist")]
fullData <- merge(all_MLE, original_data, by = "peakID")

all_genes_original = as.data.frame(table(fullData$Gene))
colnames(all_genes_original) = c("Gene.name", "Freq")
rownames(all_genes_original) = all_genes_original$Gene.name

# Mean mutations per gene
all_genes_original$MeanSub <- tapply(fullData$Nmut, fullData$Gene, mean) #Prop.Sub

# Omega
Selectome = fread("/Users/alaverre/Documents/Detecting_positive_selection/results/SQL_drosophila_omega.csv")
colnames(Selectome) = c("Ensembl.Gene.ID", "Gene.name", "omega0", "omega2", "lrt", "qvalue")
Selectome <- unique(merge(Selectome, gene_names, by = "Ensembl.Gene.ID"))

# Gene expression change
anat = "embryo"
gene.expr = fread(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/gene_expression/droso_", anat, "_stats.txt"))
gene.expr <- unique(merge(gene.expr, gene_names, by = "Ensembl.Gene.ID"))
gene.expr$abs_log2FC <- abs(gene.expr$log2FC)
fullData$meanlog2_mel <- gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$meanlog2_mel
fullData$abs_log2FC <- abs(gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$log2FC)
fullData$FDR <- -log10(gene.expr[match(fullData$Gene, gene.expr$Gene.name),]$FDR)


evol_data <- gene.expr
omega <- "abs_log2FC"
#pdf(paste0(path, "figures/drosophila_evol_", omega, "_Selectome.pdf"), width = 10, height = 6)
# Correl nb total and omega
all_genes <- merge(all_genes_original, evol_data, by = "Gene.name")
boxplot(all_genes[[omega]] ~ cut(all_genes$Freq, breaks = c(0, 10, 25, 100, 500, 1000, 5500)), notch=T, outline=F,
        xlab = "Nb Peaks", ylab = omega, main = "Nb Total Peaks vs Omega", cex=0.3)

cor <- cor.test(all_genes$Freq, all_genes[[omega]], method = "spearman")
mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)

# Ranked GO - Nb Total
geneList <- all_genes$Freq
names(geneList) <- all_genes$Gene.name
geneList <- sort(geneList, decreasing = TRUE)

gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 1, 
               maxGSSize = 500, scoreType = "pos", eps = 0)
p1 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Nb Total")

#test <- tapply(fullData$FDR_null_pos, fullData$Gene, function(x) which(p.adjust(x, method = "bonferroni")<0.5))
######################## Get Gene Lists ######################## 
genes_ID <- list()
genes_ID_background <- list()
par(mfrow=c(2,2))
for (test in c("Conclusion")){ #"Conclusion", "FDR_Conclusion", "signif", "signif.FDR"
  print(test)
  # Get Positive according to Test
  pos_genes = as.data.frame(table(fullData[fullData$FDR_Conclusion == "Directional (+)",]$Gene))
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
  
  p2 <- dotplot(go, showCategory = 20) + ggtitle("Classic GO (>= 1 Positive)")
  
  
  genes_ID[["PosGene"]] <- gene_match[gene_match$Gene.name %in% pos_genes$Gene.name,]$Gene.Stable.ID
  genes_ID_background[["PosGene"]] <- gene_match[gene_match$Gene.name %in% all_genes$Gene.name,]$Gene.Stable.ID
  
  # Ranked GO 
  geneList <- pos_genes$Freq
  names(geneList) <- pos_genes$Gene.name
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 1, 
                 maxGSSize = 500, scoreType = "pos", eps = 0)
  p3 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Nb Positive")
  
  emapplot(pairwise_termsim(gsego))
  # Ranked GO
  #go_r <- GSEA(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
  #             pvalueCutoff = 0.05, minGSSize = 1, maxGSSize = 500, pAdjustMethod = "BH")
  
  
  # Ranked GO - Ratio of positive / total
  pos_genes$all = all_genes[match(pos_genes$Gene.name, all_genes$Gene.name),]$Freq
  pos_genes_filtered <- pos_genes[pos_genes$all >= 5,]
  pos_genes_filtered$ratio = pos_genes_filtered$Freq / pos_genes_filtered$all
  boxplot(pos_genes_filtered[[omega]] ~ cut(pos_genes_filtered$ratio, breaks = c(0, 0.05, 0.1, 0.2, 1)), notch=T, 
          xlab = "Ratio Positive", ylab = omega, main=test, cex=0.3)
  
  cor <- cor.test(pos_genes_filtered$ratio, pos_genes_filtered[[omega]], method = "spearman")
  mtext(paste0("Spearman's rho = ", round(cor$estimate, 2), ", p-val = ", format.pval(cor$p.value, digits=2)), side=3, line=-2)
  
  geneList <- pos_genes_filtered$ratio
  names(geneList) <- pos_genes_filtered$Gene.name
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 1, minGSSize = 1, 
                 maxGSSize = 500, scoreType = "pos", eps = 0)
  p4 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Ratio positive/total")
  
  genes_high_ratio <- pos_genes_filtered[pos_genes_filtered$ratio >= 0.2,]$Gene.name
  
  go <- enrichGO(gene = as.character(genes_high_ratio), universe =  as.character(pos_genes_filtered$Gene.name),
                 OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 1, readable = TRUE)
  
  p5 <- dotplot(go, showCategory = 20) + ggtitle("Classic GO (Ratio >= 0.2 )")
  
  
  genes_ID[["HighRatio"]] <-  gene_match[gene_match$Gene.name %in% genes_high_ratio,]$Gene.Stable.ID
  genes_ID_background[["HighRatio"]] <- gene_match[gene_match$Gene.name %in% pos_genes_filtered$Gene.name,]$Gene.Stable.ID
  
  # write genes_ID
  #write.table(names(geneList), paste0(path, "figures/drosophila_genes_ratio_ranked_", test, ".txt"), 
  #            sep = "\t", row.names = F, col.names = F, quote = F)
  #write.table(genes_ID, paste0(path, "figures/drosophila_genes_ratio_0.4_", test, ".txt"), 
  #            sep = "\t", row.names = F, col.names = F, quote = F)
  #write.table(genes_ID_background, paste0(path, "figures/drosophila_genes_ratio_", test, ".txt"), 
  #           sep = "\t", row.names = F, col.names = F, quote = F)
  
  # Plot all figures in the same pannel
  pdf(paste0(path, "final_figures/SupFig16_GO_enrichment_", test, ".pdf"), width = 14, height = 14)
  plot <- (p1 | p2) / (p3 | p4)
  print(plot) 
  dev.off()
  
}


#dev.off()

