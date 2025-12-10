######################## Top Anat ######################## 
library(BgeeDB)

# Read Rds
pathFigure <- "/Users/alaverre/Documents/Detecting_positive_selection/final_figures/"
genes_ID <- readRDS(file = paste0(pathFigure, "/drosophila_GO_genes.rds"))
genes_ID_background <- readRDS(file = paste0(pathFigure, "/drosophila_GO_genes_background.rds"))

Enrich_classic <- list()
Enrich_weight <- list()
plots <- list()
pathResults <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"

# Create Bgee object for Drosophila
bgee <- Bgee$new(species = "7227", dataType = "rna_seq")
#stages <- c(embryo= "UBERON:0000068", larva= "UBERON:0000069", pupa= "UBERON:0000070", adult= "UBERON:0000066", all="UBERON:0000104")
stages <- c(adult= "UBERON:0000066")

for (stage in names(stages)){
  stageID = stages[[stage]]
  myTopAnatData <- loadTopAnatData(bgee, stage = stageID)
  for (Genes in c("PosGene", "HighRatio")){
    print(Genes)
    foreground <- genes_ID[[Genes]]
    background <- genes_ID_background[[Genes]]
    
    # Create gene list (1 = of interest, 0 = background)
    geneList <- as.numeric(background %in% foreground)
    names(geneList) <- background
    
    # Build TopAnat object
    myTopAnatObject <- topAnat(myTopAnatData, geneList)
    
    # Run different algorithms
    results_classic <- runTest(myTopAnatObject, algorithm = 'classic', statistic = 'fisher')
    results_weight  <- runTest(myTopAnatObject, algorithm = 'weight01', statistic = 'fisher')
    
    # Create results tables
    Enrich_classic[[stage]] <- makeTable(myTopAnatData, myTopAnatObject, results_classic, cutoff = 0.05)
    Enrich_weight[[stage]] <- makeTable(myTopAnatData, myTopAnatObject, results_weight, cutoff = 0.05)
  
  ######################## Plot TopAnat results ######################## 
  for (stats in c("Enrich_classic", "Enrich_weight")){
    # Keep significant results only
    top_terms <- get(stats)[[stage]] %>% arrange(FDR)
    
    p1 <- ggplot(top_terms, aes(x = foldEnrichment, y = reorder(organName, -log10(FDR)))) +
      geom_point(aes(size = annotated, color = -log10(FDR))) +
      scale_color_gradient(low = "lightblue", high = "red") +
      labs(x = "Fold Enrichment", y = "Anatomical Term",
           title = paste(Genes, stage, stats),
           color = "-log10(FDR)", size = "Total Genes in Term") +
      theme_minimal() +
      theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15, face = "bold")
      )
    
    print(p1)
    plots[[paste0(Genes, "_", stage, "_", stats)]] <- p1
    
    # Save plot
    ggsave(filename = paste0(pathFigure, "SupFig17_TopAnat_", Genes, "_", stats, ".png"),
           plot = p1, width = 10, height = 8, dpi=320)
    
  }
  }
}
