library(dplyr)
library(readr)
library(tidyr)
library(corrplot)

################################################################################
sp = "drosophila" # or "human"
sample = "Ni12" # or Wilson
TFs = "CTCF" # or c("CEBPA", "FOXA1", "HNF4A", "HNF6")

path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster"
pathData <- paste(path, "results/peaks_calling/NarrowPeaks", sp, sample, "bowtie2/mergedLibrary/macs2/narrowPeak/", sep="/") 

################################################################################
for (TF in TFs){
  message("######### Processing TF: ", TF, " #########")
  
  if(sp == "drosophila") {
    # ---- Single replicate workflow ----
    np_file  <- paste0(pathData, "/", TF, "_sample_peaks.narrowPeak")
    xls_file <- paste0(pathData, "/", TF, "_sample_peaks.xls")
    
    np  <- read_tsv(np_file, col_names = c(
      "chr","start","end","peak_id","score","strand",
      "signalValue","pValue","qValue","peak"))
    xls <- read_tsv(xls_file, comment = "#")
    
    # Merge and compute metrics
    final_table <- np %>%
      left_join(xls %>% select(name, pileup, fold_enrichment),
                by = c("peak_id"="name")) %>%
      mutate(
        Coverage = pileup,
        SignalStrengthScore = fold_enrichment,
        interval_id = peak_id,
        CompositeScore = log2(SignalStrengthScore + 1) * log2(Coverage + 1)
      ) %>%
      select(chr, start, end, interval_id, Coverage, SignalStrengthScore, CompositeScore) %>%
      arrange(desc(CompositeScore))
    
  }else{
    # ---- Multiple replicates workflow ----
    fc <- read_tsv(paste0(pathData, "/consensus/", TF, "/", TF, ".consensus_peaks.featureCounts.txt"), comment = "#")
    cons <- read_tsv(paste0(pathData, "/consensus/", TF, "/", TF, ".consensus_peaks.boolean.txt"))
    
    # Count replicates
    count_cols <- grep("\\.bam$", names(fc), value = TRUE)
    n_reps <- length(count_cols)
    
    # Ensure numeric columns are of numeric type
    fc_cols    <- grep("\\.fc$",  names(cons), value = TRUE)
    qval_cols <- grep("\\.qval$", names(cons), value = TRUE)
    pval_cols <- grep("\\.pval$", names(cons), value = TRUE)
    cons <- cons %>%  mutate( across(all_of(c(fc_cols, qval_cols, pval_cols)), ~ as.numeric(.)))
    
    # Normalize counts by library size (CPM)
    lib_sizes <- colSums(fc[, count_cols])
    
    fc_norm <- fc %>% mutate(
      mean_CPM = rowMeans( sweep(select(., all_of(count_cols)), 2, lib_sizes, "/") * 1e6)
    ) %>%  select(Geneid, mean_CPM)
    
    # Compute signal & consistency metrics
    signal_metrics <- cons %>%
      rowwise() %>%
      mutate(
        SignalStrengthScore = mean(c_across(all_of(fc_cols)), na.rm = TRUE),
        ConsistencyScore = {
          vals <- log2(c_across(all_of(fc_cols)))
          vals <- vals[is.finite(vals)]
          if (length(vals) <= 1) 0 else 1 / (1 + sd(vals))
        },
        SupportScore = num_samples / n_reps
      ) %>%
      ungroup()
    
    # Compute composite score
    final_table <- signal_metrics %>%
      left_join(fc_norm, by = c("interval_id" = "Geneid")) %>%
      mutate(CompositeScore = SupportScore * log2(SignalStrengthScore + 1) * ConsistencyScore
      ) %>% select(chr, start, end,interval_id, Coverage = mean_CPM, SupportScore, 
                   SignalStrengthScore, ConsistencyScore, CompositeScore
      ) %>% arrange(desc(CompositeScore))
    
    # Compute correlation matrix
    metrics <- final_table %>%  select(Coverage, SignalStrengthScore, ConsistencyScore, CompositeScore)
    cor_matrix <- cor(metrics, use = "pairwise.complete.obs")
    print(cor_matrix)
    corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.8)
    
  }
  
  # Write output
  write_tsv(final_table, paste0(pathData, "/consensus/", TF, "/", TF, ".consensus_peaks.quality_scored.tsv"))

}

################################################################################