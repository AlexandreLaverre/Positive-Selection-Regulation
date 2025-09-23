library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(ggseqlogo)
library(patchwork)
library(tidyr)

method="quantile_50bins"
all_MLE_list <- readRDS(paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/allMLE_list_", method, ".Rds"))
path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
TFS=c("HNF6", "CEBPA","HNF4A", "FOXA1")
sp="human"
sample="Wilson"
TFS = c("CEBPA")

rename_ID <- function(sp, IDs){
  if (sp == "drosophila"){
    IDs <- gsub(".*CTCF", "CTCF", IDs)
  }else{
    IDs <- gsub(".*Interval", "Interval", IDs)
  }
  return(IDs)
}

getMatrix <- function(seqs, N=250, type="SVM"){
  SVM_matrixes <- list()
  for (i in 1:N){
    seq = seqs[i,-1]
    ID = rownames(seq)
    
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
  MLE$start <- as.numeric(str_split_i(MLE$ID, ":", 2))
  MLE$length <- as.numeric(str_split_i(MLE$ID, ":", 3))-MLE$start 
  MLE$effective_length <- length[rownames(MLE)]
  MLE$deltaLength <- MLE$length-MLE$effective_length
  MLE$summits <- summits[MLE$peaks_ID,]$V2
  MLE$pos_summit <- (MLE$summits-MLE$start)+1
  MLE <- MLE[which(MLE$deltaLength==0),]
  row.names(MLE) <- rename_ID(sp, row.names(MLE))

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

  # Get matrix for each sequence
  deltaSVM_matrixes <- getMatrix(deltas_all, N=1000, type="deltaSVM")
  SVM_matrixes <- getMatrix(svm_per_base_orginal, N=1000, type="SVM")
  
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
  cor_delta_P <- c()
  cor_phylo_P <- c()
  cor_phast_P <- c()
  cor_svm_P <- c()
  lines <- list()
  summits <- c()
  max_summits <- c()
  IDS <- c()
  N = colnames(deltas)[-length(colnames(deltas))] #[1:250] or [-1]
  for (ID in N){ 
    # Ensure similar length (i.e: no indel)
    box <- boxplot(deltas[,ID]~deltas$pos, plot=F)
    deltaSVM_seq = box$stats[3,][which(!is.na(box$stats[3,]))]
    
    phyloP_seq = phyloP[which(!is.na(phyloP[,ID])),ID]
    phast_seq = phastCons[which(!is.na(phastCons[,ID])),ID]
    svm_seq = svm_per_base[which(!is.na(svm_per_base[,ID])),ID]
    
    seq_length = length(deltaSVM_seq)
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
    cor_delta[ID] <- test_delta$estimate
    cor_phylo[ID] <- test_phylo$estimate
    cor_phast[ID] <- test_phast$estimate
    
    # pearson
    test_svm <- cor.test(coverage, svm_seq, method="pearson")
    test_delta <- cor.test(coverage, -deltaSVM_seq, method="pearson")
    test_phylo <- cor.test(coverage, phyloP_seq, method="pearson")
    test_phast <- cor.test(coverage, phast_seq, method="pearson")
    cor_svm_P[ID] <- test_svm$estimate
    cor_delta_P[ID] <- test_delta$estimate
    cor_phylo_P[ID] <- test_phylo$estimate
    cor_phast_P[ID] <- test_phast$estimate
  
    dat <- data.frame(list("Normalised coverage"=coverage, "SVM"=svm_seq, 
                           "abs(deltaSVM)"=abs(deltaSVM_seq), "phastCons"=phast_seq,"phyloP"=phyloP_seq, "Position"=1:seq_length))
    lines <- c(lines, list(dat))
    IDS <- c(IDS, ID)
  }
  names(lines) <- IDS
  
  #################################################################################
  # Compare correl phastCons vs Delta
  #pdf(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/", sp, "_", TF, "_scores_correlations_along_peaks_CPM.pdf"))
  
  # Correlation Spearman
  data <- data.frame(
    name=c(rep("SVM",length(cor_svm)), rep("DeltaSVM",length(cor_delta)), rep("PhastCons",length(cor_phast)), rep("PhyloP",length(cor_phylo))),
    value=c(cor_svm, cor_delta, cor_phast, cor_phylo))

   g <- ggplot(data, aes(x=name, y=value, fill=name)) +
    geom_violin(width=1.1, alpha=0.5) +
    geom_boxplot(width=0.12, color="black", alpha=0.5) +
    scale_y_continuous(breaks = seq(-1, 1, by = 0.25)) +
    theme_minimal() +
    theme(legend.position="none", axis.title=element_text(size=14),
      plot.title = element_text(size=15), axis.text = element_text(size=12)) +
    ggtitle(paste(TF, " - Reads coverage correlations along peaks")) +
    scale_x_discrete(limits=c("SVM", "DeltaSVM", "PhastCons", "PhyloP"), ) +
    xlab("") + ylab("Spearman's correlation coefficient") +
     scale_fill_manual(values=c("navy", "red", "#7CAE00", "orange"))
   
  plot(g)
  
  # Correlation Pearson
  data <- data.frame(
    name=c(rep("SVM",length(cor_svm_P)), rep("DeltaSVM",length(cor_delta_P)), rep("PhastCons",length(cor_phast_P)), rep("PhyloP",length(cor_phylo_P))),
    value=c(cor_svm_P, cor_delta_P, cor_phast_P, cor_phylo_P))
  
  g <- ggplot(data, aes(x=name, y=value, fill=name)) +
    geom_violin(width=1.1, alpha=0.5) +
    geom_boxplot(width=0.12, color="black", alpha=0.5) +
    scale_y_continuous(breaks = seq(-1, 1, by = 0.25)) +
    theme_minimal() +
    theme(legend.position="none", axis.title=element_text(size=14),
          plot.title = element_text(size=15), axis.text = element_text(size=12)) +
    ggtitle(paste(TF, " - Reads coverage correlations along peaks")) +
    scale_x_discrete(limits=c("SVM", "DeltaSVM", "PhastCons", "PhyloP"), ) +
    xlab("") + ylab("Pearson's correlation coefficient") +
    scale_fill_manual(values=c("navy", "red", "#7CAE00", "orange"))
  plot(g)
  
  # Average centered at peak summit
  lines_summit <- list()
  size=50
  for (ID in N){
    # Ensure similar length (i.e: no indel)
    box <- boxplot(deltas[,ID]~deltas$pos, plot=F)
    deltaSVM_seq = box$stats[3,][which(!is.na(box$stats[3,]))]
    
    phyloP_seq = phyloP[which(!is.na(phyloP[,ID])),ID]
    seq_length = length(deltaSVM_seq)
    if (length(phyloP_seq) != seq_length){next}
    
    # Find summit
    coverage = unname(sum_cov[,ID])
    max_summit <- MLE[ID, "pos_summit"]
    max_summit = which(coverage==max(coverage))[1]
    end = max_summit+size
    
    # check if need to add NA
    if (max_summit<size){nbNA=size-max_summit; start = 1
    }else{nbNA=0; start = max_summit-size+1}
    
    # Sequence centred around summit
    coverage = c(rep(NA, nbNA), coverage[start:end])
    svm_seq = c(rep(NA, nbNA), svm_per_base[which(!is.na(svm_per_base[,ID])),ID][start:end])
    deltaSVM_seq =  c(rep(NA, nbNA),deltaSVM_seq[start:end])
    phast_seq =  c(rep(NA, nbNA),phastCons[which(!is.na(phastCons[,ID])),ID][start:end])
    phyloP_seq =  c(rep(NA, nbNA),phyloP_seq[start:end])
    
    dat <- data.frame(list("Normalised coverage"=coverage, "SVM"=svm_seq, 
                           "abs(deltaSVM)"=abs(deltaSVM_seq), "phastCons"=phast_seq,"phyloP"=phyloP_seq, "Position"=-(size-1):size))
    lines_summit <- c(lines_summit, list(dat))
  }
  
  Concat_lines <- do.call(rbind, lines_summit)
  Mean_lines <- as.data.frame(apply(Concat_lines, 2, function(var) tapply(var, as.factor(Concat_lines$Position), function(x) mean(x, na.rm=T))))
  melt_data <- melt(Mean_lines, id.var="Position")
  g <- ggplot(melt_data, aes(x = Position, y = value)) + geom_line(aes(color = variable)) + ylab("") +
    facet_grid(variable ~ ., scales = "free_y") + theme_bw() + theme(legend.position = "none") + 
    geom_vline(xintercept =0, linetype="dashed") +
    scale_color_manual(values = c(
      "Normalised.coverage" = "black",
      "SVM" = "orange",
      "abs.deltaSVM." = "navy", 
      "phastCons" = "red",
      "phyloP" = "#7CAE00"
    ))
  plot(g)

  pdf(paste0("/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/final_figures/Fig1_A_", sp, "_", TF, "_scores_CPM.pdf"))
  
  # Plot 20 peaks example
  for (ID in N){
    if (ncol(SVM_matrixes[[ID]])>150){next}
    if (cor_svm[ID]<0.5){next}
    message(ID)
    custom_labels <- c("Normalised.coverage" = "Coverage", "SVM" = "SVM", "abs.deltaSVM." = "ΔSVM", "phastCons" = "phastCons", "phyloP" = "phyloP")
    ID = names(deltaSVM_matrixes[ID])
    #if (ID != "Interval_43502"){next}
    melt_data <- melt(lines[[ID]], id.var="Position")
    #melt_data <- melt_data[which(melt_data$variable %in% c("Normalised.coverage", "phastCons", "phyloP")),]
    line_plot <- ggplot(melt_data, aes(x = Position, y = value)) + geom_line(aes(color = variable)) + ylab("") +
      facet_grid(variable ~ ., scales = "free_y", switch = "y", labeller = labeller(variable = custom_labels)) + theme_bw() + 
      theme(legend.position = "none", strip.placement = "outside", strip.text.y = element_text(size = 11), strip.background = element_blank()) +
      geom_vline(xintercept = max_summits[ID], linetype="dashed") +
      scale_color_manual(values = c(
        "Normalised.coverage" = "black",
        "SVM" = "orange", 
        "abs.deltaSVM." = "navy", 
        "phastCons" = "red",
        "phyloP" = "#7CAE00"
      ))
    #plot(line_plot)

    SVM_seqlogo <- ggseqlogo(SVM_matrixes[[ID]], method = "custom") + theme_minimal() +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank()) + 
      labs(y = "SVM", x = "") + ggtitle(ID) + scale_x_continuous(breaks = seq(0, size-1, by = 25))
    #plot(SVM_seqlogo)
    
    deltaSVM_seqlogo <- ggseqlogo(deltaSVM_matrixes[[ID]], method = "custom") + theme_minimal() +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank()) + 
      labs(y = "deltaSVM", x = "") + ggtitle("") + scale_x_continuous(breaks = seq(0, size-1, by = 25))
    tmp_plot <- SVM_seqlogo / deltaSVM_seqlogo 

    final_plot <- (SVM_seqlogo / deltaSVM_seqlogo / line_plot) +
      plot_layout(nrow = 3, heights = c(2, 2, 10)) &
      theme(plot.margin = margin(1, 1, 1, 1))
    
    print(final_plot)
    
  }
  
  dev.off()
}

