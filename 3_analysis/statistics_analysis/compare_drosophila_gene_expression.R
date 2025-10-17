library(BgeeDB)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

species <- c("Drosophila_melanogaster", "Drosophila_simulans")

orthologs <- read.csv("/Users/alaverre/Documents/Detecting_positive_selection/docs/orthologs_7227-7240.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(orthologs) <- c("mel_gene", "sim_gene")
orthologs <- orthologs[, c("mel_gene", "sim_gene")]
ortho <- orthologs %>% group_by(mel_gene) %>% filter(n() == 1) %>%        # keep only mel_genes that appear once
  ungroup() %>% group_by(sim_gene) %>% filter(n() == 1) %>%               # also ensure each sim_gene appears once
  ungroup()

avg.exp <- list()
exp.pivot <- list()
for (sp in species) {
  print(paste("Species:", sp))
  bgee <- Bgee$new(species = sp, dataType = "rna_seq")

  for (anat in c("embryo", "adult")){
    # Embryo
    if (anat == "embryo"){
      anatEntityId = "UBERON:0000922"
      stageID = "UBERON:0000068"
    }else{
      #Whole adult organism : days 5-8 of fully formed stage 
      anatEntityId = "UBERON:0007023"
      stageID = ifelse(sp == "Drosophila_melanogaster", "FBdv:00007080", "DsimDv:0000007")
    }
    
    data <- getSampleProcessedData(bgee, anatEntityId = anatEntityId, stageId = stageID)
    message(sp, ", ", anat, "; Number of samples retrieved: ", length(unique(data$Library.ID)))
    expr <- as.data.frame(formatData(bgee, data, stats = "tpm")@assayData[["exprs"]])
    expr$gene_id <- row.names(expr)
    avg.exp[[sp]][[anat]] <- expr %>% mutate(mean_TPM = rowMeans(across(-gene_id), na.rm = TRUE)) %>%  select(gene_id, mean_TPM)
    exp.pivot[[sp]][[anat]] <-  expr %>% pivot_longer(-gene_id, names_to = "sample", values_to = "TPM") %>% mutate(species = sp)
    
  }
}

################################################################################
for (anat in c("embryo", "adult")){
  # Average TPM per gene for each species
  merged <- ortho %>%
    inner_join(avg.exp[["Drosophila_melanogaster"]][[anat]], by = c("mel_gene" = "gene_id")) %>%
    inner_join(avg.exp[["Drosophila_simulans"]][[anat]], by = c("sim_gene" = "gene_id"),
               suffix = c("_mel", "_sim"))
  
  merged <- merged %>% filter(mean_TPM_mel > 1 | mean_TPM_sim > 1)
  
  merged <- merged %>% mutate(
    logTPM_mel = log2(mean_TPM_mel + 1),
    logTPM_sim = log2(mean_TPM_sim + 1),
    logTPM_mel_centered = logTPM_mel - median(logTPM_mel, na.rm = TRUE),
    logTPM_sim_centered = logTPM_sim - median(logTPM_sim, na.rm = TRUE),
    delta_log2 = abs(logTPM_mel_centered - logTPM_sim_centered)
  )
  
  
  summary_stats <- merged %>%
    summarise(
      pearson_r = cor(logTPM_mel_centered, logTPM_sim_centered, method = "pearson"),
      spearman_r = cor(logTPM_mel_centered, logTPM_sim_centered, method = "spearman"),
      median_delta = median(delta_log2),
      mean_delta = mean(delta_log2)
    )
  
  print(summary_stats)
  # Scatter plot of expression (mel vs sim)
  p1 <- ggplot(merged, aes(x = logTPM_mel_centered, y = logTPM_sim_centered)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = anat,
      x = "Mean expression in D. melanogaster",
      y = "Mean expression in D. simulans"
    ) +
    theme_minimal(base_size = 14)
  
  ################################################################################
  # Volcano plot with all replicates
  
  expr_long <- bind_rows(exp.pivot[["Drosophila_melanogaster"]][[anat]], exp.pivot[["Drosophila_simulans"]][[anat]])
  expr_long <- expr_long %>% select(gene_id, species, sample, TPM)
  expr_long_filtered <- expr_long %>%  filter(
    (species == "Drosophila_melanogaster" & gene_id %in% ortho$mel_gene) |
      (species == "Drosophila_simulans" & gene_id %in% ortho$sim_gene)
  )
  expr_long_filtered <- expr_long_filtered %>% filter(!is.na(TPM) )
  expr_long_filtered <- expr_long_filtered %>% filter(TPM > 0.5) 
  expr_long_filtered <- expr_long_filtered %>% mutate(logTPM = log2(TPM + 1))
  
  # Normalise logTPM by centering on median per species
  expr_long_filtered <- expr_long_filtered %>%
    group_by(species) %>%
    mutate(logTPM_centered = logTPM - median(logTPM, na.rm = TRUE)) %>%
    ungroup()
  
  expr_long_filtered <- expr_long_filtered %>%
    mutate(Ensembl.Gene.ID = case_when(
      species == "Drosophila_melanogaster" ~ gene_id,
      species == "Drosophila_simulans" ~ ortho$mel_gene[match(gene_id, ortho$sim_gene)]
    )
    ) %>%
    filter(!is.na(Ensembl.Gene.ID))  # keep only 1:1 orthologs
  
  stats <- expr_long_filtered %>%
    group_by(Ensembl.Gene.ID) %>%
    summarise(meanlog2_mel = mean(logTPM_centered[species == "Drosophila_melanogaster"]),
              meanlog2_sim = mean(logTPM_centered[species == "Drosophila_simulans"]),
              log2FC = meanlog2_mel - meanlog2_sim,
              pval = tryCatch(t.test(logTPM_centered ~ species)$p.value, error = function(e) NA)
    ) %>%
    ungroup() %>%
    mutate(
      negLogP = -log10(pval),
      significant = pval < 0.01 & abs(log2FC) > 1
    ) 
  
  stats$FDR <- p.adjust(stats$pval, method = "fdr")
  stats$FDR_signif <- stats$FDR < 0.1 & abs(stats$log2FC) > 1
  
  p2 <- ggplot(stats, aes(x = log2FC, y = -log10(FDR), color = FDR_signif)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("gray70", "red"), name = "Significance", guide="none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgray") +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "darkgray") +
    labs(
      title = "",
      x = "log2 fold-change",
      y = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 14)
  
  combined <- p1 + p2
  
  # save stats table
  write.table(stats, file = paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/gene_expression/droso_", anat, "_stats.txt"), sep = "\t", quote = F, row.names = F)
  ggsave(paste0("/Users/alaverre/Documents/Detecting_positive_selection/final_figures/SupFig_drosophilas_expression", anat, ".pdf"), width = 12, height = 6)
  
}


