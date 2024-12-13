library(progress)
args = commandArgs(trailingOnly=TRUE)

path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
sp = ifelse(is.na(args[1]), "mouse", args[1]) 
TF = ifelse(is.na(args[1]), "CEBPA", args[2]) 
maxSub=150

################################################################################
### DATAS
deltas <- list()
for (data in c("focal")){ #"sister"
  # all possible deltas
  all.deltas.file = paste0(path, "positive_selection/all_deltas/", sp, "/", TF, "/", "all_possible_deltaSVM.txt") #data
  
  all_IDs = paste(rep(paste0("pos", 1:1000), each=3), c("sub1", "sub2", "sub3"), sep = ":") # because max seq length=1000
  all_col = c("seq_name", all_IDs)
  all.deltas <- read.table(all.deltas.file, h=F, sep="\t", quote="", fill=T, col.names = all_col)
  row.names(all.deltas) <- all.deltas$seq_name
  message("Nb sequences ", data, ": ", nrow(all.deltas))
  deltas[[paste0("all.deltas.", data)]] <- all.deltas[,2:3001]
  
  # observed deltas
  obs.deltas.file = paste0(path, "positive_selection/all_deltas/", sp, "/", TF, "/", "observed_deltaSVM.txt") #data
  
  obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
  obs <- read.table(obs.deltas.file, h=F, sep="\t", quote="", fill=T, col.names = obs_col)
  row.names(obs) <- obs$seq_name
  deltas[[paste0("obs.deltas.", data)]] <- obs[rownames(all.deltas), 5:105]
}

################################################################################
seq.focal <- rownames(deltas[["obs.deltas.focal"]])
#seq.sister <- rownames(deltas[["obs.deltas.sister"]])
#common.seq <- intersect(seq.focal, seq.sister)

####### Simulations 
permut = 1000
seq.nb = 10000 #length(seq.focal)  #length(common.seq)

empty.df = data.frame(matrix(nrow = seq.nb, ncol = permut), row.names = seq.focal[1:seq.nb])
deltas.simul <- list("delta_drift"=empty.df, "delta_purif_abs_sister"=empty.df,
                     "delta_purif_abs_focal"=empty.df, "delta_pos"=empty.df, "delta_neg"=empty.df )

pval.higher <- obs[seq.focal[1:seq.nb],2:4]

# Create a progress bar
pb <- progress_bar$new(total = seq.nb, format = "[:bar] :percent :elapsed")

##### Run permutations
for (seq in seq.focal[1:seq.nb]){
  all.muts = deltas[["all.deltas.focal"]][seq,]
  all.muts = all.muts[!is.na(all.muts)]
  nb.sub = pval.higher[seq, "NbSub"]
  delta.obs = pval.higher[seq, "deltaSVM"]
  
  #sister.muts = deltas[["obs.deltas.sister"]][seq,]
  #sister.muts = sister.muts[!is.na(sister.muts)]
  
  focal.muts = deltas[["obs.deltas.focal"]][seq,]
  focal.muts = focal.muts[!is.na(focal.muts)]
  
  #abs.sis = max(abs(sister.muts))
  abs.foc = max(abs(focal.muts))
  
  #mut.abs.sis = all.muts[which(all.muts>-abs.sis & all.muts<abs.sis)]
  mut.abs.foc = all.muts[which(all.muts>-abs.foc & all.muts<abs.foc)]
  
  mut.pos = all.muts[which(all.muts>0)]
  # mut.neg = all.muts[which(all.muts<0)]
    
  for (i in 1:permut){
    deltas.simul[["delta_drift"]][seq,i] <- sum(sample(all.muts, nb.sub))
    deltas.simul[["delta_purif_abs_focal"]][seq,i] <- sum(sample(mut.abs.foc, nb.sub))
    deltas.simul[["delta_pos"]][seq,i] <- sum(sample(mut.pos, nb.sub))
    #deltas.simul[["delta_purif_abs_sister"]][seq,i] <- sum(sample(mut.abs.sis, nb.sub))
    #deltas.simul[["delta_neg"]][seq,i] <- sum(sample(mut.neg, nb.sub))
  }
  
  pval.higher[seq,"drift"] <- sum(deltas.simul[["delta_drift"]][seq,] > delta.obs) / permut
  pval.higher[seq,"purif_abs_focal"] <- sum(deltas.simul[["delta_purif_abs_focal"]][seq,] > delta.obs) / permut
  pval.higher[seq,"pos"] <- sum(deltas.simul[["delta_pos"]][seq,] > delta.obs) / permut
  #pval.higher[seq,"purif_abs_sister"] <- sum(deltas.simul[["delta_purif_abs_sister"]][seq,] > ddelta.obs) / permut
  #pval.higher[seq,"neg"] <- sum(deltas.simul[["delta_neg"]][seq,] > delta.obs) / permut
  
  pb$tick()
}


pval.higher$type0.05 <- ifelse(pval.higher$drift<0.05 & pval.higher$purif_abs_focal<0.05, "positive", 
                           ifelse(pval.higher$drift<0.05 & pval.higher$purif_abs_focal>0.05, "only-non-drift", "drift"))
pval.higher$type0.01 <- ifelse(pval.higher$drift<0.01 & pval.higher$purif_abs_focal<0.01, "positive", 
                               ifelse(pval.higher$drift<0.01 & pval.higher$purif_abs_focal>0.01, "only-non-drift", "drift"))

save(list=c("pval.higher", "deltas.simul", "deltas"), file=paste0(path, "tested_deltas_", sp, "_", TF, ".Rdata"))


################################################################################
## TPM PLOTS
sister <- as.matrix(deltas[["obs.deltas.sister"]])
sister <- sister[!is.na(sister)]

focal <- as.matrix(deltas[["obs.deltas.focal"]])
focal <- focal[!is.na(focal)]

all.sister <- as.matrix(deltas[["all.deltas.sister"]])
all.sister <- all.sister[!is.na(all.sister)]

all.focal <- as.matrix(deltas[["all.deltas.focal"]])
all.focal <- all.focal[!is.na(all.focal)]

boxplot(focal, sister, all.focal, all.sister, notch=T, outline=F,
        ylab="deltaSVM per substitution",
        names=c("obs.focal", "obs.sister", "all.focal", "all.sister"))
abline(h=0)

################################################################################
hist(obs$SVM, breaks=100)
obs$SVM.class <- cut(obs$SVM, breaks=c(min(obs$SVM), -400, -200, -150, -75, max(obs$SVM)))
boxplot(abs(obs$deltaSVM)~obs$SVM.class, notch=T, outline=F)
boxplot(abs(obs$deltaSVM)~as.factor(obs$NbSub), outline=F, notch=T, xlim=c(0,21), ylim=c(0,30))
boxplot(obs$deltaSVM~as.factor(obs$NbSub), outline=F, notch=T, xlim=c(0,21), ylim=c(-30,20))
