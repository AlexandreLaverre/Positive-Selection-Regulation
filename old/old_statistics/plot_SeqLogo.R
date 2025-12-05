library(ggseqlogo)
library(tidyr)
library(ggplot2)
library(stringr)

path = "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/positive_selection/"
pathFigures <- "/Users/alaverre/Documents/Detecting_positive_selection/results/figures/"
TF="CEBPA"
sp="human"
sample="Wilson"

################################################################################
# 1 - Load data 
# deltas
file_all_deltas = paste0(path, "/NarrowPeaks/human/Wilson/", TF, "/deltas/ancestral_all_possible_deltaSVM.txt")
deltas_all <- read.table(file_all_deltas, h=T, sep="\t", quote="", fill=T)
deltas_Interval <- gsub(".*:Interval", "Interval", deltas_all$ID)
rownames(deltas_all) <- deltas_Interval
deltas_all <- as.data.frame(deltas_all[,-1])
deltas <- as.data.frame(t(deltas_all))
deltas$pos <- sub("\\..*", "", rownames(deltas))
deltas$pos <- factor(deltas$pos, levels=str_sort(unique(deltas$pos), numeric = TRUE))

# SVM
all_SVM = paste0(path, "/NarrowPeaks/", sp, "/", sample, "/", TF, "/SVM_per_base.txt")
all_SVM <- read.table(all_SVM, h=T, sep="\t", quote="", fill=T)
rownames(all_SVM) <-  gsub(".*:Interval", "Interval", all_SVM$ID)
all_SVM <- as.data.frame(all_SVM[,-1])

################################################################################
# 2 - Make sequence logos for SVM and deltaSVM
getMatrix <- function(seq){
  # Convert wide to long format using pivot_longer
  long_df <- seq %>% pivot_longer(cols = everything(),
                                  names_to = c("Position", "Nucleotide"), names_sep = "\\.", values_to = "Score")
  
  # Convert back into matrix format using pivot_wider
  matrix_df <- long_df %>% pivot_wider(names_from = Position, values_from = Score)
  nuc = matrix_df[1]
  matrix_df <- as.matrix(matrix_df[-1])
  rownames(matrix_df) <- as.vector(nuc)[[1]]
  
  # Remove position without any information
  matrix_df <- matrix_df[, colSums(!is.na(matrix_df)) > 0]
  seq_len = length(colnames(matrix_df))
  
  return(matrix_df)
}

SVM_matrixes <- list()
deltaSVM_matrixes <- list()
for (i in 1:250){
  ID = rownames(deltas_all[i,])
  
  # Use deltaSVM to retrieve sequence and ensure no indel
  delta_seq <- deltas_all[ID,]
  actual_seq <- colnames(delta_seq[which(is.na(delta_seq))])[1:4000]
  delta_noNA <- boxplot(deltas[,ID]~deltas$pos, plot=F)
  seq_length_delta = length(delta_noNA$stats[3,][which(!is.na(delta_noNA$stats[3,]))])
  names <- colnames(delta_seq[which(is.na(delta_seq))])[1:1000]
  
  # Extract SVM 
  seq_svm <- all_SVM[ID,]
  names(seq_svm) <- names
  seq_length = length(seq_svm[which(!is.na(seq_svm))])
  
  if (seq_length != seq_length_delta){next}
  
  SVM_matrixes[[ID]] <- getMatrix(seq_svm)
  deltaSVM_matrixes[[ID]] <- getMatrix(delta_seq)
}

################################################################################
# 3 - Plot sequence logos
for (i in 1:250){
  if (ncol(SVM_matrixes[[i]])>150){next}
  print(i)
  ID = names(SVM_matrixes[i])
  seqlogo <- ggseqlogo(SVM_matrixes[[ID]], method = "custom") + theme_minimal() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text.x = element_blank()) +
    labs(y = plot, x = "Positions") + scale_x_continuous(breaks = seq(0, 150-1, by = 25)) + 
    ggtitle(ID)
  
  plot(seqlogo)
}

################################################################################