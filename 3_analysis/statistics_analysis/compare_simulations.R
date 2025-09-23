# Compare simulations 

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/BroadPeaks/human/Wilson/CEBPA/Tests/"

type="positive"
# All simulations files
permut.file <- paste0(path, "PosSelTest_deltaSVM_10000permutations_simulation_by_params_independent_SVM_without_backMut_noBin_", type, ".txt")
permut_epistasis.file <- paste0(path, "PosSelTest_deltaSVM_10000permutations_simulation_by_params_epistasis_without_backMut_noBin_", type, ".txt")
MLE_old.file <- paste0(path, "old/new_run_low_w.csv") 
MLE_new.file <- paste0(path, "old/new_run_low_w_onlymut.csv") 
hist_100bins.file <- paste0(path, "MLE_summary_simulated_last_", type, "_hist_100bins_threshold_0.01.csv")
hist_50bins.file <- paste0(path, "MLE_summary_simulated_last_", type, "_hist_50bins_threshold_0.01.csv")
quant_100bins.file <- paste0(path, "MLE_summary_simulated_last_", type, "_quantile_100bins_threshold_0.01.csv")
quant_50bins.file <- paste0(path, "MLE_summary_simulated_last_", type, "_quantile_50bins_threshold_0.01.csv")


# Get data in a single list
simul <- c("permut", "permut_epistasis", "MLE_old", "MLE_new", "hist_100bins", "hist_50bins", "quant_100bins", "quant_50bins")

par(mfrow=c(2,2))
simulations <- list()
for (data in simul){
  if (grepl("permut", data)){
    simulations[[data]] <- read.csv(get(paste0(data, ".file")), h=T, row.names = 1, sep="\t")
    order=row.names(simulations[[data]])
    simulations[[data]]$Conclusion <- ifelse(simulations[[data]]$pval.high<0.01, "Directional (+)",
                                             ifelse(simulations[[data]]$pval.high>0.99, "Directional (-)", "Neutral model"))
  }else{
    simulations[[data]] <- read.csv(get(paste0(data, ".file")), h=T, row.names = 1)
    simulations[[data]] <- simulations[[data]][order,]
  }
  
  IDs = rownames(simulations$permut)
  simulations[[data]] = simulations[[data]][IDs,]
  plot(table(simulations[[data]]$Conclusion), main=data, ylab="Nb simul")
}

par(mfrow=c(1,1))
boxplot(simulations[["quant_50bins"]]$Nmut~simulations[["quant_50bins"]]$Conclusion, main="Simul Directional - Quantile 50bins",
        xlab="Conclusion", ylab="Nb Substitutions")

boxplot(simulations[["quant_50bins"]]$MeanObs~simulations[["quant_50bins"]]$Conclusion, main="Simul Directional - Quantile 50bins",
        xlab="Conclusion", ylab="Mean Obs")

test <- simulations[["quant_100bins"]]
test$Dir <- ifelse(test$Conclusion=="Neutral model", 0, ifelse(test$Conclusion=="Disruptive", 0, 1))
plot(prop.table(table(test$Dir,test$Nmut), margin = 2)[2,]~seq(2,10,1), main="Simul Directional - Hist 100bins",
     xlab="Nb Sub", ylab="Proportion of Directional", type="b")


# All deltas 
Alldelta.file <- paste0(pathResults, "/deltas/focal_all_possible_deltaSVM.txt")
Alldeltas <- read.table(Alldelta.file, h=T, sep="\t", quote="", fill=T, row.names = 1)
#Alldeltas <- Alldeltas[order,]

# Observed deltas
delta.file <- paste0(path, "/deltas/simul_deltas_positive_observed_deltaSVM.txt")
obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
deltas <- read.table(delta.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
row.names(deltas) <- deltas$seq_name
deltas <- deltas[rownames(Alldeltas),]
row.names(rand) <- rand$ID
row.names(pos) <- pos$ID

IDs = row.names(test[which(test$Conclusion=="Neutral model" & test$Nmut >4),])[1:16]
IDs = row.names(rand[which(rand$Conclusion=="Directional (-)" & rand$Nmut >4),])[1:4]
IDs = row.names(pos[which(pos$Conclusion=="Disruptive"),])[1:4]
par(mfrow=c(2,2))
for (ID in IDs){
  #ID="chr2:31837530:31837664"
  Alpha=pos[ID,]$p_value_null_pos
  hist(as.numeric(Alldeltas[ID,]), breaks=200,  main=paste("pval Null/Pos=",signif(Alpha, 2)), xlab="deltaSVM")
  abline(v=as.numeric(deltas[ID,5:14]), col="red")
}


