library(dplyr)
library(VennDiagram)
library(ggplot2)
library(ggvenn)

options("stringsAsFactors"=FALSE)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/peaks_overlap/"
TF <- "CEBPA"
species <- c("Mmusculus", "Mspretus", "Mcaroli")

################################################################################
peaks <- list()
lifted_peaks <- list()
for (sp in species){
  BED <- read.table(paste0(path, "peaks/", sp, "_", TF, ".peaks.bed"))
  colnames(BED) <- c("chr", "start", "end", "complete_ID")
  row.names(BED) <- paste0(BED$chr, ":", BED$start, ":", BED$end)
  BED[,sp] <- row.names(BED)
  
  for (target in setdiff(species, sp)){
    # Correspondence lifted
    lifted <- read.table(paste0(path, "/lifted_peaks/", sp, "_", TF, "_to_", target, "_0.9.bed"), header=F)
    lifted[[target]] <- paste0(lifted$V1, ":", lifted$V2, ":", lifted$V3)
    row.names(lifted) <- lifted$V4
    lifted_peaks[[sp]][[target]] <- lifted
    
    # Overlapping peaks
    overlap <- read.table(paste0(path, "/", TF, "_", sp, "_lifted_overlap_", target, ".txt"), header=T)
    BED[overlap$ID, target] <- overlap$overlap_ID
  }
  
  peaks[[sp]] <- BED[,species]
}

combined <- bind_rows(peaks)
row.names(combined) <- NULL
data <- unique(combined)
row.names(data) <- 1:nrow(data)

################################################################################
# Create a boolean data.frame for Venn diagramm
data_bool <- replace(data, is.na(data), FALSE)
data_bool <- replace(data_bool, !is.na(data), TRUE)

data_bool$Mmusculus <- as.logical(data_bool$Mmusculus)
data_bool$Mcaroli <- as.logical(data_bool$Mcaroli)
data_bool$Mspretus <- as.logical(data_bool$Mspretus)


data_bool$Nb <- 1:nrow(data_bool)
sets_list <- list(Mmusculus = data_bool[which(data_bool$Mmusculus), "Nb"],
                  Mcaroli = data_bool[which(data_bool$Mcaroli), "Nb"],
                  Mspretus = data_bool[which(data_bool$Mspretus), "Nb"])

#venn.plot <- venn.diagram(sets_list, filename = NULL)
#grid.draw(venn.plot)
#fit <- euler(data_bool)
#plot(fit)

ggvenn(sets_list)

################################################################################
# Write output
for (sp in species){
  conserved_row <- row.names(data_bool[which(data_bool$Mmusculus == TRUE & data_bool$Mcaroli == TRUE & data_bool$Mspretus == TRUE),])
  conserved <- data[conserved_row, sp]
  
  other_sp1 = setdiff(species, sp)[1]
  other_sp2 = setdiff(species, sp)[2]
  gain_row <- row.names(data_bool[which(data_bool[[sp]] == TRUE & 
                                          data_bool[[other_sp1]] == FALSE &
                                          data_bool[[other_sp2]] == FALSE),])
  
  loss_row <- row.names(data_bool[which(data_bool[[sp]] == FALSE &
                                          data_bool[[other_sp1]] == TRUE &
                                          data_bool[[other_sp2]] == TRUE),])
  
  # Get the correspondence IDs from other_sp to sp
  loss_IDs <- data[loss_row, other_sp1]
  loss_in_sp <- lifted_peaks[[other_sp1]][[sp]][loss_IDs, sp]
  
  write.table(data[conserved_row, sp], file = paste0(path, "/peaks/", sp, "_conserved_peaks.txt"), quote=F, col.names = F, row.names = F)
  write.table(data[gain_row, sp], file = paste0(path, "/peaks/", sp, "_gain_peaks.txt"), quote=F, col.names = F, row.names = F)
  write.table(loss_in_sp, file = paste0(path, "/peaks/", sp, "_loss_peaks.txt"), quote=F, col.names = F, row.names = F)
}


