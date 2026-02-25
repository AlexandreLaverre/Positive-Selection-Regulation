library("matrixStats")
library("data.table")
set.seed(123)

path <- "/Users/alaverre/Documents/Detecting_positive_selection/"

sp="human"
Permut_file = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/CTCF_binding/CTCF/Tests/PosSelTest_deltaSVM_10000permutations_last.txt")
RegEvol_file = paste0(path, "cluster/results/positive_selection/NarrowPeaks/", sp, "/CTCF_binding/CTCF/Tests/MLE_summary_exact.csv")

treshold = 0.1
minMut = 2

dataFilter<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  if (ncol(deltaSVM) >= 7){
    deltaSVM <- deltaSVM[,c("ID", "deltaSVM", "NbSub", "pval.high")]
    nsplit=3;minus=0
  }else{nsplit=5; minus=1}
  
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  deltaSVM <- deltaSVM[which(deltaSVM$varNumb >= minMut),]
  
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=nsplit, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-minus
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:ncol(deltaSVM))])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$ID <- paste(ID.bed[,1],ID.bed[,2],ID.bed[,3], sep="_")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  
  deltaSVM$FDR<-p.adjust(deltaSVM$pValue, method="fdr")
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  
  rownames(deltaSVM) <- deltaSVM$ID
  return(deltaSVM)
}

# Data
Permut <- dataFilter(fread(Permut_file))

RegEvol <- read.csv(RegEvol_file, row.names = 1)
RegEvol <- deltaSVM_ML[which(deltaSVM_ML$Nmut >= minMut),]
RegEvol$FDR <- p.adjust(deltaSVM_ML$p_value_purif_pos, method="fdr")

list_deltaSVM <- list(Permut, RegEvol)
names(list_deltaSVM) <- c("Permut", "RegEvol")

### Peaks annotation
organNumb<-fread(paste0(path, "/cluster/data/human_ENCODE_CTCF_merged_annotated.bed"))
## from the fourth column, if the value > 0, indicating this binding site is expressed in this tissue
tissue <- list()
tissue[["human"]] <- c("adrenal gland","B cell","esophagus muscularis mucosa","retinal pigment epithelial cell",
                       "omental fat pad","gastrocnemius medialis","astrocyte of the cerebellum","astrocyte of the spinal cord",
                       "brain microvascular endothelial cell","choroid plexus epithelial cell","heart left ventricle","liver","upper lobe of left lung",
                       "neural cell (in vitro differentiated)","ovary", "pancreas","foreskin fibroblast","peyer patch","prostate gland","lower leg skin",
                       "spleen","stomach","testis","thoracic aorta","thyroid gland","tibial nerve","transverse colon","uterus","vagina")

tissue[["mouse"]] <- c("Bone marrow","Cerebellum","Cortical plate","Heart","Kidney","Liver",
                       "Lung","Olfactory bulb","Small intestine","Testis","Thymus")

# Define organs to systems and its colors
col <- c("#E41A1C", "#FF7F00", "#66C2A5", "#999999", "#FFFF33", "#F781BF", "#A65628", "#984EA3", "#377EB8", "black", "navy")
names(col) <- c("Nervous", "Male reproductive","Immune", "Endocrine", "Integumentary",
                "Respiratory", "Cardiovascular", "Digestive", "Female reproductive", "Skeletomuscular", "Excretory")

systems <- list()
systems[["human"]] <- c("Endocrine", "Immune", "Digestive", "Nervous", "Integumentary", "Skeletomuscular",
                        "Nervous", "Nervous", "Nervous", "Nervous", "Cardiovascular", "Endocrine", 
                        "Respiratory", "Nervous", "Female reproductive", "Endocrine", "Male reproductive", "Immune",
                        "Male reproductive", "Integumentary", "Immune", "Digestive", "Male reproductive", "Cardiovascular",
                        "Endocrine", "Nervous", "Digestive", "Female reproductive", "Female reproductive")

systems[["mouse"]] <- c("Immune", "Nervous", "Nervous", "Cardiovascular", 
                        "Excretory", "Endocrine", "Respiratory", "Nervous", 
                        "Digestive", "Male reproductive", "Immune")

#systems <- paste(systems, "system")
organs <- colnames(organNumb)[c(4:Ncol)]
names(systems[[sp]]) <- organs

colnames(organNumb)<-c("chr","start","end", tissue[[sp]])
Ncol=ncol(organNumb)
organNumb$orgNumb<-apply(organNumb[,4:Ncol], 1, function(x) length(x[x>0]))


deltaSVMorgNumb <- merge(deltaSVM, organNumb,by=c("chr","start","end"))
saveRDS(deltaSVMorgNumb, paste0(path, "human_CTCF_ENCODE_RegEvol.Rds"))

resample_n <- 5000
n_resample <- 10000

# Initialize lists to store results
SumDeltaLL_resampled <- list()
MedianDeltaLL_resampled <- list()
posNumb_resampled <- list()

# Compute Proportion and SUM of deltaLL for each tissue with resampling
for (i in tissue[[sp]]) {
  message(paste("Tissue:", i))
  tissueBasedData <- deltaSVMorgNumb[which(deltaSVMorgNumb[, i] > 0), ]
  
  # Remove tissues with less than resample_n peaks
  if (nrow(tissueBasedData) < resample_n){
    message(paste("Removing tissue:", i, "with", nrow(tissueBasedData), "peaks"))
    tissue[[sp]] <- tissue[[sp]][tissue[[sp]] != i]
    next
  }
  
  # RESAMPLING: draw 10,000 subsets of 5,000 peaks
  SumDeltaLL_resampled[[i]] <- numeric(n_resample)
  for (r in 1:n_resample){
    resample_data <- tissueBasedData[sample(nrow(tissueBasedData), resample_n), ]
    resample_data$deltaLL <- 2 * (resample_data$LL_pos - resample_data$LL_neutral)
    posNumb_resampled[[i]][r] <- sum(resample_data$FDR < threshold)
    SumDeltaLL_resampled[[i]][r] <- sum(resample_data$deltaLL^(1/4))  # 4th-root as in Daub et al.
    MedianDeltaLL_resampled[[i]][r] <- median(resample_data$deltaLL^(1/4))  # 4th-root as in Daub et al.
  }
  
}


# Create a data frame for plotting
plotDF <- data.frame(
  tissue = tissue[[sp]],
  system = systems[[sp]][tissue[[sp]]],
  
  SumDeltaLL = sapply(SumDeltaLL_resampled, median),
  sd_SumDeltaLL = sapply(SumDeltaLL_resampled, sd),
  
  MedianDeltaLL = sapply(MedianDeltaLL_resampled, median),
  sd_MedianDeltaLL = sapply(MedianDeltaLL_resampled, sd),
  
  posNumb = sapply(posNumb_resampled, median),
  sd_posNumb = sapply(posNumb_resampled, sd)
)

stat="SumDeltaLL"

pdf(paste0(path, "final_figures/", sp, "_CTCF_SUMSTAT.pdf"), width=9, height=8)
par(mfrow=c(1,1))
par(mar=c(10.5, 5, 2, 9.5) + 0.1)

plotDF <- plotDF[order(plotDF[[stat]]), ]
Norgan <- nrow(plotDF)

plot(1:Norgan, plotDF[[stat]], pch=16, cex=2,
     ylim=c(min(plotDF[[stat]] - plotDF[[paste0("sd_", stat)]]),
            max(plotDF[[stat]] + plotDF[[paste0("sd_", stat)]])),
     col=col[plotDF$system], xaxt="n",
     xlab="", ylab="SUMSTAT", main="", cex.lab=1.2, cex.main=1.5, cex.axis=1.2)

# Add error bars for resampled SD
arrows(1:Norgan, plotDF[[stat]] - plotDF[[paste0("sd_", stat)]],
       1:Norgan, plotDF[[stat]] + plotDF[[paste0("sd_", stat)]],
       angle=90, code=3, length=0.05, col=col[plotDF$system])

# Add vertical grid
abline(v=1:Norgan, col="grey", lty=4)

# Add tissue labels
axis(side = 1, at = 1:Norgan, labels=F)
text(1:Norgan, par("usr")[3]-3, srt = 45, cex=1, adj=1, labels = plotDF$tissue, xpd=TRUE)

# Add legend
legend("topright", legend=names(col), title = "System",title.font = 2, col = col, cex=1, pch=16, pt.cex=2, bty="n", inset=c(-0.32,0), xpd=TRUE)

dev.off()

# Test Nervous system vs other
for (i in unique(systems[["human"]])){
  nervous_tissues <- names(systems[[sp]])[systems[[sp]] == i]
  other_tissues <- setdiff(tissue[[sp]], nervous_tissues)
  
  # Combine resampled SUM scores
  nervous_sum <- unlist(SumDeltaLL_resampled[nervous_tissues])
  other_sum <- unlist(SumDeltaLL_resampled[other_tissues])
  
  # Wilcoxon rank-sum test (one-sided: Nervous > Others)
  wilcox_test <- wilcox.test(nervous_sum, other_sum, alternative = "greater")
  cliff_test <- cliff.delta(nervous_sum, other_sum)
  
  # Print results
  print(paste("Comparing", i, "vs Others: Wilcoxon p-value =", wilcox_test$p.value, 
              "Cliff's delta =", cliff_test$estimate, "magnitude=", cliff_test$magnitude))
  
}



