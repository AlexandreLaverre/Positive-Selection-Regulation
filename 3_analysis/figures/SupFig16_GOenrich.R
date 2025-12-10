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
pathFigure <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
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

all_genes = as.data.frame(table(fullData$Gene))
colnames(all_genes) = c("Gene.name", "Freq")
rownames(all_genes) = all_genes$Gene.name

################################################################################
test="FDR_Conclusion"
genes_ID <- list()
genes_ID_background <- list()

# Ranked GO - Nb Total
geneList <- all_genes$Freq
names(geneList) <- all_genes$Gene.name
geneList <- sort(geneList, decreasing = TRUE)

gsego <- gseGO(geneList = geneList, OrgDb = "org.Dm.eg.db", keyType = "SYMBOL", ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, minGSSize = 1, 
               maxGSSize = 500, scoreType = "pos", eps = 0)
p1 <- dotplot(gsego, showCategory = 20, split=".sign") + ggtitle("Ranked List - Nb Total")


# Get Positive according to Test
pos_genes = as.data.frame(table(fullData[fullData[[test]] == "Directional (+)",]$Gene))
colnames(pos_genes) = c("Gene.name", "Freq")
rownames(pos_genes) = pos_genes$Gene.name

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

# Ranked GO - Ratio of positive / total
pos_genes$all = all_genes[match(pos_genes$Gene.name, all_genes$Gene.name),]$Freq
pos_genes_filtered <- pos_genes[pos_genes$all >= 5,]
pos_genes_filtered$ratio = pos_genes_filtered$Freq / pos_genes_filtered$all

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
               pAdjustMethod = "BH", qvalueCutoff = 1, readable = TRUE)

p5 <- dotplot(go, showCategory = 20) + ggtitle("Classic GO (Ratio >= 0.1)")


genes_ID[["HighRatio"]] <-  gene_match[gene_match$Gene.name %in% genes_high_ratio,]$Gene.Stable.ID
genes_ID_background[["HighRatio"]] <- gene_match[gene_match$Gene.name %in% pos_genes_filtered$Gene.name,]$Gene.Stable.ID

# Plot all figures in the same pannel
plot <- (p1 | p2) / (p3 | p5)
print(plot) 

# Save plot
ggsave(plot = plot, filename = paste0(pathFigure, "/SupFig16_GO_enrich.pdf"), width = 12, height = 14, dpi = 320)
ggsave(plot = plot, filename = paste0(pathFigure, "/SupFig16_GO_enrich.png"), width = 12, height = 14, dpi = 320)

# Save Gene Lists as Rds
saveRDS(genes_ID, file = paste0(pathFigure, "/drosophila_GO_genes.rds"))
saveRDS(genes_ID_background, file = paste0(pathFigure, "/drosophila_GO_genes_background.rds"))

