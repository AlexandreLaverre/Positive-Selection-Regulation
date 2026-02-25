library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(progress)

path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
all_MLE_list <- readRDS(paste0("/Users/alaverre/Documents/Detecting_positive_selection/ChIP_peaks_RegEvol_exact_ranked.Rds"))

sp="human"
sample="Wilson"
TFS=c("CEBPA")
rename_ID <- function(sp, IDs){
  if (sp == "drosophila"){
    IDs <- gsub(".*CTCF", "CTCF", IDs)
  }else{
    IDs <- gsub(".*Interval", "Interval", IDs)
  }
  return(IDs)
}

getMatrix <- function(seqs, IDs, type="SVM"){
  SVM_matrixes <- list()
  for (ID in IDs){
    seq = seqs[ID,-1]
    
    # Use deltaSVM to retrieve sequence
    if (type=="SVM"){
      delta_seq <- deltas_all[ID,-1]
      names <- colnames(delta_seq[which(is.na(delta_seq))])[1:1000]
      names(seq) <- names
    }
    
    # Convert wide to long format using pivot_longer
    long_df <- seq %>% pivot_longer(cols = everything(), names_to = c("Position", "Nucleotide"), names_sep = "\\.", values_to = "Score")
    
    # Convert back into matrix format using pivot_wider
    matrix_df <- long_df %>% pivot_wider(names_from = Position, values_from = Score)
    nuc = matrix_df[1]
    matrix_df <- as.matrix(matrix_df[-1])
    rownames(matrix_df) <- as.vector(nuc)[[1]]
    
    # Remove position without any information
    matrix_df <- matrix_df[, colSums(!is.na(matrix_df)) > 0]
    
    SVM_matrixes[[ID]] <- matrix_df
  }
  
  return(SVM_matrixes)
  
}

getScore <- function(score, sp, sample, TF){
  file = paste0(path, "../", score, "/NarrowPeaks/", sp, "/", sample, "/", TF, "/score_per_base.txt")
  score_base <- read.table(file, h=F, sep="\t", quote="", fill=T, col.names = paste0("pos", 0:3000), row.names = 1)
  # remove row with more than 1500 non-NA values
  score_base <- score_base[rowSums(!is.na(score_base)) <= 1500, 1:1500]
  rownames(score_base) <- rename_ID(sp, rownames(score_base))
  score_base <- score_base[row.names(MLE),]
  score_base <- as.data.frame(t(score_base))
  return(score_base)
} 

#all_data <- list()
for (TF in TFS){
  message(TF)
  ################################################################################
  # Peak summits
  file_summits = paste0(path, "../peaks_calling/NarrowPeaks/", sp, "/", sample, "/", TF, ".consensus_summits_UCSC_names.bed")
  summits <- read.table(file_summits, h=F, sep="\t", quote="", fill=T)
  summits$V4 <- gsub(".*:\\d+:\\d+_", "", summits$V4)
  rownames(summits) <- summits$V4
  
  # Sequence details
  sequences <- readDNAStringSet(paste0(path, "NarrowPeaks/", sp,  "/", sample, "/", TF, "/sequences/filtered_focal_sequences.fa"), format = "fasta")
  length <- width(sequences)
  names(length) <- gsub("_.*:\\d+:\\d+:", "_", names(sequences))
  MLE <- all_MLE_list[[paste0(sp, "_", TF)]]
  
  if (sp == "drosophila"){rownames(MLE) <- paste0(MLE$peaks_ID)
  }else{rownames(MLE) <- paste0(MLE$ID, "_", MLE$peaks_ID)}
  
  MLE$start <- as.numeric(str_split_i(MLE$ID, ":", 2))
  MLE$length <- as.numeric(str_split_i(MLE$ID, ":", 3))-MLE$start 
  
  sep <- ifelse(sp == "drosophila", ":", "_")
  MLE$effective_length <- length[paste0(MLE$ID, sep, MLE$peaks_ID)]
  MLE$deltaLength <- MLE$length-MLE$effective_length
  MLE$summits <- summits[MLE$peaks_ID,]$V2
  MLE$pos_summit <- (MLE$summits-MLE$start)+1
  MLE <- MLE[which(MLE$deltaLength<5),]
  MLE$meanSVM <- MLE$SVM / MLE$effective_length
  row.names(MLE) <- rename_ID(sp, row.names(MLE))
  
  # Peaks quality
  pathData <- paste(path, "../peaks_calling/NarrowPeaks", sp, sample, "bowtie2/mergedLibrary/macs2/narrowPeak/consensus", sep="/") 
  peaks_quality <- read_tsv(paste0(pathData, "/", TF, "/", TF, ".consensus_peaks.quality_scored.tsv"), show_col_types = FALSE)
  peaks_quality <- as.data.frame(peaks_quality)
  rownames(peaks_quality) <- peaks_quality$interval_id
  peaks_quality <- peaks_quality[row.names(MLE),]
  peaks_quality$mean_coverage <- peaks_quality$Coverage / MLE$effective_length
  
  # SVM
  svm_per_base = paste0(path, "/NarrowPeaks/", sp, "/", sample, "/", TF, "/SVM_per_base.txt")
  svm_per_base <- read.table(svm_per_base, h=T, sep="\t", quote="", fill=T)
  rownames(svm_per_base) <-   rename_ID(sp, svm_per_base$ID)
  svm_per_base_orginal <- svm_per_base[row.names(MLE),]
  svm_per_base <- as.data.frame(t(svm_per_base[,-1]))
  
  # deltas
  file_all_deltas = paste0(path, "/NarrowPeaks/", sp,  "/", sample, "/", TF, "/deltas/ancestral_all_possible_deltaSVM.txt")
  deltas_all <- read.table(file_all_deltas, h=T, sep="\t", quote="", fill=T)
  deltas_Interval <- rename_ID(sp, deltas_all$ID)
  rownames(deltas_all) <- deltas_Interval
  deltas_all <- deltas_all[row.names(MLE),]
  
  deltas <- as.data.frame(t(deltas_all[,-1]))
  deltas$pos <- sub("\\..*", "", rownames(deltas))
  deltas$pos <- factor(deltas$pos, levels=str_sort(unique(deltas$pos), numeric = TRUE))
  
  # PhyloP and PhastCons
  phyloP <- getScore("phyloP", sp, sample, TF)
  phastCons <-  getScore("phastCons", sp, sample, TF)
  
  # Reads Coverage
  pathCoverage = paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/peaks_calling/NarrowPeaks/", sp, "/", sample, "/bowtie2/mergedLibrary/deepTools/coverage/")
  coverage_files <- list.files(paste0(pathCoverage,"/", TF), pattern = "\\matrix_CPM.gz$")
  coverage <- list()
  for (file in coverage_files){
    cov <- read.table(gzfile(paste(pathCoverage, TF, file, sep="/")), skip=1)
    row.names(cov) <- cov$V4
    cov <- cov[row.names(MLE),]
    coverage[[file]] <- t(cov[,7:ncol(cov)])
  } 
  
  if (sp == "drosophila"){
    mean_cov <-  coverage[[file]]
    sum_cov <-  coverage[[file]]
  }else{
    mean_cov <- (coverage[[1]]+coverage[[2]])/2
    sum_cov <- coverage[[1]]+coverage[[2]] 
  }
  
  ################################################################################
  # Compute correlation and plots
  cor_delta <- c()
  cor_phylo <- c()
  cor_phast <- c()
  cor_svm <- c()
  cor_svm_pval <- c()
  cor_delta_P <- c()
  cor_phylo_P <- c()
  cor_phast_P <- c()
  cor_svm_P <- c()
  cor_svm_P_pval <- c()
  lines <- list()
  summits <- c()
  max_summits <- c()
  IDS <- c()
  diff = c()
  N = colnames(deltas)[-length(colnames(deltas))] #[1:250] or [-1]
  
  pb <- progress_bar$new(
    total = length(N),
    format = "Processing [:bar] :percent | ETA: :eta"
  )
  
  for (ID in N) {
    # Ensure similar length (i.e: no indel)
    box <- boxplot(deltas[,ID]~deltas$pos, plot=F)
    deltaSVM_seq = box$stats[3,][which(!is.na(box$stats[3,]))]
    
    phyloP_seq = phyloP[which(!is.na(phyloP[,ID])),ID]
    phast_seq = phastCons[which(!is.na(phastCons[,ID])),ID]
    svm_seq = svm_per_base[which(!is.na(svm_per_base[,ID])),ID]
    
    seq_length = length(deltaSVM_seq)
    diff = c(diff, length(phyloP_seq)-seq_length)
    if (length(phyloP_seq) != seq_length){next}
    
    #message(colnames(phyloP)[i])
    coverage = sum_cov[,ID][1:seq_length]
    
    summits[ID] <- MLE[ID, "pos_summit"]
    max_summits[ID] <- which(coverage==max(coverage))[1]
    
    # Correlations
    test_svm <- cor.test(coverage, svm_seq, method="spearman")
    test_delta <- cor.test(coverage, -deltaSVM_seq, method="spearman")
    test_phylo <- cor.test(coverage, phyloP_seq, method="spearman")
    test_phast <- cor.test(coverage, phast_seq, method="spearman")
    cor_svm[ID] <- test_svm$estimate
    cor_svm_pval[ID] <- test_svm$p.value
    cor_delta[ID] <- test_delta$estimate
    cor_phylo[ID] <- test_phylo$estimate
    cor_phast[ID] <- test_phast$estimate
    
    # pearson
    test_svm <- cor.test(coverage, svm_seq, method="pearson")
    test_delta <- cor.test(coverage, -deltaSVM_seq, method="pearson")
    test_phylo <- cor.test(coverage, phyloP_seq, method="pearson")
    test_phast <- cor.test(coverage, phast_seq, method="pearson")
    cor_svm_P[ID] <- test_svm$estimate
    cor_svm_P_pval[ID] <- test_svm$p.value
    cor_delta_P[ID] <- test_delta$estimate
    cor_phylo_P[ID] <- test_phylo$estimate
    cor_phast_P[ID] <- test_phast$estimate
    
    dat <- data.frame(list("Normalised coverage"=coverage, "SVM"=svm_seq, 
                           "abs(deltaSVM)"=abs(deltaSVM_seq), "phastCons"=phast_seq,"phyloP"= phyloP_seq, "Position"=1:seq_length))
    lines <- c(lines, list(dat))
    IDS <- c(IDS, ID)
    pb$tick()
  }
  
  names(lines) <- IDS
  peaks_quality$cor_svm_P <- cor_svm_P[row.names(peaks_quality)]
  peaks_quality$cor_svm_P_pval <- cor_svm_P_pval[row.names(peaks_quality)]
  peaks_quality$cor_svm <- cor_svm[row.names(peaks_quality)]
  
  peaks_quality$cor_delta_P <- cor_delta_P[row.names(peaks_quality)]
  peaks_quality$cor_phylo_P <- cor_phylo_P[row.names(peaks_quality)]
  peaks_quality$cor_phast_P <- cor_phast_P[row.names(peaks_quality)]
  
  MLE$class_length <- cut(MLE$effective_length, breaks=quantile(MLE$effective_length, probs=seq(0,1, by=0.1), na.rm=TRUE))
  
  peaks_quality$class_coverage <- cut(peaks_quality$Coverage, breaks=quantile(peaks_quality$Coverage, probs=seq(0,1, by=0.1), na.rm=TRUE), labels=seq(1,10), , include.lowest = TRUE)
  peaks_quality$class_SignalStrengthScore <- cut(peaks_quality$SignalStrengthScore, breaks=quantile(peaks_quality$SignalStrengthScore, probs=seq(0,1, by=0.1), na.rm=TRUE), labels=seq(1,10), , include.lowest = TRUE)
  all_data[[paste0(sp, "_", TF)]] <- peaks_quality
}


# Save RDS file for all_data
saveRDS(all_data, "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/Correlation_SVM_Coverage_Quality.Rds")
