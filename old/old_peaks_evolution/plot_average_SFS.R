path <- "/Users/alaverre/Documents/Detecting_positive_selection/cluster/results/polymorphism_analyses/NarrowPeaks/drosophila/modERN/"

# Get all SFS 
samples <- gsub(paste0(path, "/"), "", list.dirs(path, recursive = FALSE))
all_SFS <- list()
for (sample in samples){
  print(sample)
  file = paste(path, sample, "SFS.csv", sep="/")
  if (!file.exists(file)){next}
  SFS <- read.csv(file, sep="\t", row.names=1)
  all_SFS[[sample]] <- SFS
}

# Plot mean SFS
mean_SFS <- log2(Reduce("+", all_SFS) / length(all_SFS))
rownames(mean_SFS) <- paste0("S=", 0:(nrow(mean_SFS)-1))

colors <- colorRampPalette(c("blue",  "red"))(nrow(mean_SFS))
matplot(t(mean_SFS), type = "l", lty = 1, col = colors, las=1,
        xlab = "Derived Allele Count", ylab = "Proportion of mutations (log10)", main = "Mean SFS")
legend("topright", legend = rownames(mean_SFS), col = colors, lty = 1, cex = 0.6)

size <- c()
nrow <- c()
for (i in names(all_SFS)){
  nrow <- c(nrow, nrow(all_SFS[[i]]))
  size <- c(size, length(all_SFS[[i]]))
}


