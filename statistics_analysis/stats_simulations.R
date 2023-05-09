library(RColorBrewer)
library(data.table)

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species <- c("dog", "human", "macaca", "mouse", "rat", "chicken", "cat","rabbit")
TFs <- c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

#pdf(paste0(path, "/figures/human_simulations.pdf"), width=8)
##############################################################################
# Distribution nb substitutions
NbSub <- list()
col <- pal[1:length(species)]
names(col) <- species

boxcol <- c()
for (sp in species){
  for (TF in TFs){
    file = paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt")
    if (file.exists(file)){
      test <- read.table(file,h=T)
      #hist(test$NbSub, xlim=c(0, 30), breaks=200, main=paste(sp,TF), xlab = "Nb substitution")
      NbSub[[paste(sp,TF, sep="_")]] <- test$NbSub
      boxcol <- c(boxcol, col[sp])
    }
  }
}

boxplot(NbSub, notch=T, outline=F, col=boxcol, cex=0.1,  xlab="Samples", ylab="Substitutions per peak")
legend("topleft", fill=col, legend=species, ncol=2, bty='n', )

##############################################################################
# Simulation: detection of positive set according to number of substitution

Nb_simul_Sub <- c(seq(2,10,1), seq(20, 100, 10))

Simulations <- list()
for (nbSimul in Nb_simul_Sub){
  simul = paste0(path, "positive_selection/human/simulation/Positive_Selection_test/PosSelTest_deltaSVM_", nbSimul, "_mutations.txt")
  Simulations[[nbSimul]] <- read.table(simul, h=T)
}


par(mfrow=(c(2,2)))
par(mai=c(0.6,0.8,0.4,0.3), mgp=c(2.2,0.8,0))
nbSub <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$NbSub)
boxplot(nbSub, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="Nb substitutions")
abline(h=c(10, 20, 40, 60, 80, 100), lty=3)

DeltaSVM <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$deltaSVM)
boxplot(DeltaSVM, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="deltaSVM")

pval <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$pval.high)
boxplot(pval, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="pval high")

propSignif <- lapply(Nb_simul_Sub, function(x) length(which(Simulations[[x]]$pval.high<=0.05))/nrow(Simulations[[x]]))
boxplot(propSignif, names=Nb_simul_Sub, xlab="Nb mutations", ylab="Proportion pval<0.05")



hist(Simulations[[3]]$deltaSVM, breaks=50, xlab="deltaSVM", main="3 mutational steps", xlim=c(-25, 40))
hist(Simulations[[3]]$med.deltaSVM.simul, breaks=50, xlab="random deltaSVM", main="3 mutational steps", xlim=c(-25, 40))

hist(Simulations[[10]]$deltaSVM, breaks=100, xlab="deltaSVM", main="10 mutational steps", xlim=c(-25, 40))
hist(Simulations[[10]]$med.deltaSVM.simul, breaks=50, xlab="random deltaSVM", main="10 mutational steps", xlim=c(-25, 40))
#dev.off()

##############################################################################
# Simulation oriented delta
Nb_simul_Sub <- c(2, 3, 4, 10, 20, 50)

par(mfrow=(c(2,3)))
par(mai=c(0.6,0.8,0.4,0.3), mgp=c(2.2,0.8,0))
Simulations_selection <- list()
for (nbSimul in Nb_simul_Sub){
  simul = read.table(paste0(path, "positive_selection/human/simulation/Positive_Selection_test/PosSelTest_deltaSVM_", nbSimul, "_mutations_selection.txt"), h=T)
  plot(simul$deltaSVM~simul$pval.high, cex=0.1, xlab="Pval.high", ylab="deltaSVM", main=paste(nbSimul, "mutational steps"))
  abline(v=0.05, col="red")
  Simulations_selection[[nbSimul]] <- simul
}

delta_class <- lapply(Nb_simul_Sub, function(x) cut(Simulations_selection[[x]]$deltaSVM, breaks=c(-50, -1, 1, 50), include.lowest=T))

par(mfrow=(c(2,2)))
par(mai=c(0.6,0.8,0.4,0.3), mgp=c(2.2,0.8,0))
nbSub <- lapply(Nb_simul_Sub, function(x) Simulations_selection[[x]]$NbSub)
boxplot(nbSub, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="Nb substitutions")
abline(h=c(10, 20, 40, 60, 80, 100), lty=3)

DeltaSVM <- lapply(Nb_simul_Sub, function(x) Simulations_selection[[x]]$deltaSVM)
boxplot(DeltaSVM, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="deltaSVM")

pval <- lapply(Nb_simul_Sub, function(x) Simulations_selection[[x]]$pval.high)
boxplot(pval, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations", ylab="pval high")

propSignif <- lapply(Nb_simul_Sub, function(x) length(which(Simulations_selection[[x]]$pval.high<=0.05))/nrow(Simulations_selection[[x]]))
boxplot(propSignif, names=Nb_simul_Sub, xlab="Nb mutations", ylab="Proportion pval<0.05", ylim=c(0.05,0.1))
