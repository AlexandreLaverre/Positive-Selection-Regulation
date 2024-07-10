library(Biostrings)
library(stringr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(reshape2)

all_MLE_list <- readRDS("/Users/alaverre/Documents/Detecting_positive_selection/results/allMLE_list.Rds")
path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
TFS=c("FOXA1", "HNF4A", "HNF6", "CEBPA")

for (TF in TFS){
  message(TF)
    
  ################################################################################
  # Peak summits
  file_summits = paste0(path, "../peaks_calling/NarrowPeaks/human/Wilson/", TF, ".consensus_summits_UCSC_names.bed")
  summits <- read.table(file_summits, h=F, sep="\t", quote="", fill=T)
  summits$V4 <- gsub(".*:\\d+:\\d+_", "", summits$V4)
  rownames(summits) <- summits$V4
  
  # Sequence details
  sequences <- readDNAStringSet(paste0(path, "NarrowPeaks/human/Wilson/", TF, "/sequences/filtered_focal_sequences.fa"), format = "fasta")
  length <- width(sequences)
  names(length) <- gsub("_.*:\\d+:\\d+:", "_", names(sequences))
  human <- all_MLE_list[[paste0("human_", TF)]]
  human$start <- as.numeric(str_split_i(human$ID, ":", 2))
  human$length <- as.numeric(str_split_i(human$ID, ":", 3))-human$start 
  human$effective_length <- length[rownames(human)]
  human$deltaLength <- human$length-human$effective_length
  human$summits <- summits[human$peaks_ID,]$V2
  human$pos_summit <- (human$summits-human$start)+1
  human <- human[which(human$deltaLength==0),]
  row.names(human) <- gsub(".*_Interval", "Interval", row.names(human))
  
  # SVM
  svm_per_base = paste0(path, "/NarrowPeaks/human/Wilson/", TF, "/SVM_per_base.txt")
  svm_per_base <- read.table(svm_per_base, h=T, sep="\t", quote="", fill=T)
  rownames(svm_per_base) <-  gsub(".*:Interval", "Interval", svm_per_base$ID)
  svm_per_base <- svm_per_base[row.names(human),]
  svm_per_base <- as.data.frame(t(svm_per_base[,-1]))
  
  # deltas
  file_all_deltas = paste0(path, "/NarrowPeaks/human/Wilson/", TF, "/deltas/ancestral_all_possible_deltaSVM.txt")
  deltas_all <- read.table(file_all_deltas, h=T, sep="\t", quote="", fill=T)
  deltas_Interval <- gsub(".*:Interval", "Interval", deltas_all$ID)
  rownames(deltas_all) <- deltas_Interval
  deltas_all <- deltas_all[row.names(human),]
  deltas <- as.data.frame(t(deltas_all[,-1]))
  deltas$pos <- sub("\\..*", "", rownames(deltas))
  deltas$pos <- factor(deltas$pos, levels=str_sort(unique(deltas$pos), numeric = TRUE))
  
  # PhyloP and PhastCons
  getScore <- function(score){
    file = paste0(path, "../", score, "/NarrowPeaks/human/Wilson/", TF, "/score_per_base.txt")
    score_base <- read.table(file, h=F, sep="\t", quote="", fill=T, col.names = paste0("pos", 0:1500), row.names = 1)
    rownames(score_base) <- gsub(":Inter", "_Inter", rownames(score_base))
    rownames(score_base) <- gsub(".*_Interval", "Interval", rownames(score_base))
    score_base <- score_base[row.names(human),]
    score_base <- as.data.frame(t(score_base))
    return(score_base)
  } 
  
  phyloP <- getScore("phyloP")
  phastCons <-  getScore("phastCons")
  
  # Reads Coverage
  pathCoverage = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/peaks_calling/NarrowPeaks/human/Wilson/bowtie2/mergedLibrary/deepTools/coverage/"
  coverage_files <- list.files(paste0(pathCoverage,"/", TF), pattern = "\\matrix_readCount.gz$")
  coverage <- list()
  for (file in coverage_files){
    cov <- read.table(gzfile(paste(pathCoverage, TF, file, sep="/")), skip=1)
    row.names(cov) <- cov$V4
    cov <- cov[row.names(human),]
    coverage[[file]] <- t(cov[,7:ncol(cov)])
  } 
  
  mean_cov <- (coverage[[1]]+coverage[[2]])/2
  sum_cov <- coverage[[1]]+coverage[[2]]
  
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
  for (i in seq(1,ncol(deltas)-1)){ 
    # Ensure similar length (i.e: no indel)
    box <- boxplot(deltas[,i]~deltas$pos, plot=F)
    deltaSVM_seq = box$stats[3,][which(!is.na(box$stats[3,]))]
    
    phyloP_seq = phyloP[which(!is.na(phyloP[,i])),i]
    phast_seq = phastCons[which(!is.na(phastCons[,i])),i]
    svm_seq = svm_per_base[which(!is.na(svm_per_base[,i])),i]
    
    seq_length = length(deltaSVM_seq)
    if (length(phyloP_seq) != seq_length){next}
    
    #message(colnames(phyloP)[i])
    coverage = sum_cov[,i][1:seq_length]
    
    ID = colnames(deltas)[i]
    summits <- c(summits, human[ID, "pos_summit"])
    max_summits <- c(max_summits, which(coverage==max(coverage))[1])
    
    # Correlations
    test_svm <- cor.test(coverage, svm_seq, method="spearman")
    test_delta <- cor.test(coverage, -deltaSVM_seq, method="spearman")
    test_phylo <- cor.test(coverage, phyloP_seq, method="spearman")
    test_phast <- cor.test(coverage, phast_seq, method="spearman")
    cor_svm <- c(cor_svm, test_svm$estimate)
    cor_delta <- c(cor_delta, test_delta$estimate)
    cor_phylo <- c(cor_phylo, test_phylo$estimate)
    cor_phast <- c(cor_phast, test_phast$estimate)
    
    # pearson
    test_svm <- cor.test(coverage, svm_seq, method="pearson")
    test_delta <- cor.test(coverage, -deltaSVM_seq, method="pearson")
    test_phylo <- cor.test(coverage, phyloP_seq, method="pearson")
    test_phast <- cor.test(coverage, phast_seq, method="pearson")
    cor_svm_P <- c(cor_svm_P, test_svm$estimate)
    cor_delta_P <- c(cor_delta_P, test_delta$estimate)
    cor_phylo_P <- c(cor_phylo_P, test_phylo$estimate)
    cor_phast_P <- c(cor_phast_P, test_phast$estimate)
  
    dat <- data.frame(list("Normalised coverage"=coverage, "SVM"=svm_seq, 
                           "abs(deltaSVM)"=abs(deltaSVM_seq), "phastCons"=phast_seq,"phyloP"=phyloP_seq, "Position"=1:seq_length))
    lines <- c(lines, list(dat))
  }
  
  
  #################################################################################
  # Compare correl phastCons vs Delta
  pdf(paste0("/Users/alaverre/Documents/Detecting_positive_selection/results/figures/", TF, "_scores_correlations_along_peaks_readCount.pdf"))
  
  # Correlation Spearman
  data <- data.frame(
    name=c(rep("SVM",length(cor_svm)), rep("DeltaSVM",length(cor_phylo)), rep("PhastCons",length(cor_phylo)), rep("PhyloP",length(cor_phylo))),
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
    xlab("") + ylab("Spearman's correlation coefficient")
   
  plot(g)
  
  # Correlation Pearson
  data <- data.frame(
    name=c(rep("SVM",length(cor_svm)), rep("DeltaSVM",length(cor_phylo)), rep("PhastCons",length(cor_phylo)), rep("PhyloP",length(cor_phylo))),
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
    xlab("") + ylab("Pearson's correlation coefficient")
  
  plot(g)
  
  # Average centered at peak summit
  lines_summit <- list()
  size=200
  for (i in seq(1,ncol(deltas)-1)){
    # Ensure similar length (i.e: no indel)
    box <- boxplot(deltas[,i]~deltas$pos, plot=F)
    deltaSVM_seq = box$stats[3,][which(!is.na(box$stats[3,]))]
    
    phyloP_seq = phyloP[which(!is.na(phyloP[,i])),i]
    seq_length = length(deltaSVM_seq)
    if (length(phyloP_seq) != seq_length){next}
    
    # Find summit
    coverage = unname(sum_cov[,i][1:seq_length])
    ID = colnames(deltas)[i]
    max_summit <- human[ID, "pos_summit"]
    #max_summit = which(coverage==max(coverage))[1]
    end = max_summit+size
    
    # check if need to add NA
    if (max_summit<size){nbNA=size-max_summit; start = 1
    }else{nbNA=0; start = max_summit-size+1}
    
    # Sequence centred around summit
    coverage = c(rep(NA, nbNA), coverage[start:end])
    svm_seq = c(rep(NA, nbNA), svm_per_base[which(!is.na(svm_per_base[,i])),i][start:end])
    deltaSVM_seq =  c(rep(NA, nbNA),deltaSVM_seq[start:end])
    phast_seq =  c(rep(NA, nbNA),phastCons[which(!is.na(phastCons[,i])),i][start:end])
    phyloP_seq =  c(rep(NA, nbNA),phyloP_seq[start:end])
    
    dat <- data.frame(list("Normalised coverage"=coverage, "SVM"=svm_seq, 
                           "abs(deltaSVM)"=abs(deltaSVM_seq), "phastCons"=phast_seq,"phyloP"=phyloP_seq, "Position"=-199:size))
    lines_summit <- c(lines_summit, list(dat))
  }
  
  Concat_lines <- do.call(rbind, lines_summit)
  Mean_lines <- as.data.frame(apply(Concat_lines, 2, function(var) tapply(var, as.factor(Concat_lines$Position), function(x) mean(x, na.rm=T))))
  melt_data <- melt(Mean_lines, id.var="Position")
  g <- ggplot(melt_data, aes(x = Position, y = value)) + geom_line(aes(color = variable)) + ylab("") +
    facet_grid(variable ~ ., scales = "free_y") + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept =0, linetype="dashed")
  plot(g)
  
  # Plot 20 peaks example
  for (i in 10:20){
    melt_data <- melt(lines[[i]], id.var="Position")
    g <- ggplot(melt_data, aes(x = Position, y = value)) + geom_line(aes(color = variable)) + ylab("") +
      facet_grid(variable ~ ., scales = "free_y") + theme_bw() + theme(legend.position = "none") + geom_vline(xintercept = max_summits[i], linetype="dashed")
    plot(g)
  }
  
  dev.off()
}

