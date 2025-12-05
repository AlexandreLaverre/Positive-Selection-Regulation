######################## Top Anat ######################## 
library(BgeeDB)

#genes_ID <- 
#genes_ID_background <- 


Enrich_classic <- list()
Enrich_elim <- list()
Enrich_weight <- list()
plots <- list()
pathResults <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"

# Create Bgee object for Drosophila
bgee <- Bgee$new(species = "7227", dataType = "rna_seq")
stages <- c(embryo= "UBERON:0000068", larva= "UBERON:0000069", pupa= "UBERON:0000070", adult= "UBERON:0000066", all="UBERON:0000104")
stages <- c(all= "UBERON:0000066")

for (stage in names(stages)){
  stageID = stages[[stage]]
  myTopAnatData <- loadTopAnatData(bgee, stage = stageID)
  for (test in c("HighRatio")){
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
    table_classic <- makeTable(myTopAnatData, myTopAnatObject, results_classic, cutoff = 0.1)
    table_elim    <- makeTable(myTopAnatData, myTopAnatObject, results_elim, cutoff = 0.1)
    table_weight  <- makeTable(myTopAnatData, myTopAnatObject, results_weight, cutoff = 0.1)
    
    # save
    if (length(table_classic) > 1){
      Enrich_classic[[stage]] <- table_classic
    }
    if (length(table_elim) >1){
      Enrich_elim[[stage]] <- table_elim
      Enrich_weight[[stage]] <- table_weight # [,c(2,6,8)]
    }
    
  }
  
  ######################## Plot TopAnat results ######################## 
  for (stats in c("Enrich_classic", "Enrich_elim", "Enrich_weight")){
    # Keep significant results only
    top_terms <- get(stats)[[stage]] %>%
      filter(FDR < 0.05) %>%
      arrange(FDR) 
    
    p1 <- ggplot(top_terms, aes(x = foldEnrichment, y = reorder(organName, -log10(pValue)))) +
      geom_point(aes(size = annotated, color = -log10(FDR))) +
      scale_color_gradient(low = "lightblue", high = "red") +
      labs(x = "Fold Enrichment", y = "Anatomical Term",
           title = paste(stage, stats),
           color = "-log10(FDR)", size = "Total Genes in Term") +
      theme_minimal() +
      theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15, face = "bold")
      )
    
    print(p1)
    plots[[paste(stage, "_", stats)]] <- p1
  }
}

# Save results
save(Enrich_classic, Enrich_elim, Enrich_weight, plots, file = paste0(pathResults, "drosophila_TopAnat_results.RData"))