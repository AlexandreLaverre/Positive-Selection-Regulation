library(RColorBrewer)
library(data.table)

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 

path <- "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species <- c("dog", "macaca", "mouse", "cat", "human", "rat", "rabbit", "chicken")
TFs <- c("CTCF", "CEBPA", "FOXA1", "HNF4A", "HNF6")

pdf(paste0(path, "/figures/human_simulations.pdf"), width=8)

##############################################################################
# Distribution nb substitutions
NbSub <- list()
NbSub_total <- list()
col_sp <- pal[1:length(species)]
names(col_sp) <- species

col_TF <- pal[1:length(TFs)]
names(col_TF) <- TFs

boxcol_sp <- c()
boxcol_TF <- c()
for (sp in species){
  for (TF in TFs){
    file = paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt")
    if (file.exists(file)){
      test <- read.table(file,h=T)
      #hist(test$NbSub, xlim=c(0, 30), breaks=200, main=paste(sp,TF), xlab = "Nb substitution")
      NbSub[[paste(sp,TF, sep="_")]] <- test$NbSub
      boxcol_sp <- c(boxcol_sp, col_sp[sp])
      boxcol_TF <- c(boxcol_TF, col_TF[TF])
      
      NbSub_total[[sp]] <- c(NbSub_total[[sp]], test$NbSub)
    }
  }
}

a <- boxplot(NbSub_total, notch=F, outline=F, col=col_sp[species],
             las=2, ylab="Substitutions per peak", cex.lab=1.2, cex.axis=1.2)

a <- boxplot(NbSub, notch=T, outline=F, col=boxcol_sp, cex=0.1, xaxt='n', las=2,
        xlab="Samples", ylab="Substitutions per peak")

#legend("topleft", fill=col_sp, legend=species, ncol=2, bty='n', )

axis(1, at = c(3, 8, 13, 18, 23, 26, 27.5, 29), labels = NA, cex.axis = 1, srt=45)
text(x = c(3, 8, 13, 18, 23, 26, 27.5, 29), y = par("usr")[3] - 5,
     labels = c("dog", "human", "macaca", "mouse", "rat", "chicken", "cat", "rabbit"),
     srt = 45, adj = 1, xpd = TRUE, cex=1)


boxplot(NbSub, notch=F, outline=F, col=boxcol_TF, cex=0.1, xaxt='n', las=2,
        xlab="Samples", ylab="Substitutions per peak")
legend("topleft", fill=col_TF, legend=TFs, ncol=2, bty='n', )

##############################################################################
# Simulation: detection of positive set according to number of substitution

Nb_simul_Sub <- c(seq(2,10,1), seq(20, 100, 10))
#Nb_simul_Sub <- c(2,4,8,20)

Simulations <- list()
for (nbSimul in Nb_simul_Sub){
  simul_file = paste0(path, "positive_selection/human/simulation_mutational_steps/Positive_Selection_test/PosSelTest_deltaSVM_", nbSimul, "_mutations.txt")
  simul <- read.table(simul_file, h=T)
  simul$FDR <- p.adjust(simul$pval.high, method="fdr")
  Simulations[[nbSimul]] <- simul
}

pdf(paste0(path, "/figures/human_100_simulations_no_selection.pdf"), width=8)
par(mfrow=(c(2,2)))
par(mai=c(0.6,0.8,0.4,0.3), mgp=c(2.2,0.8,0))

DeltaSVM <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$deltaSVM)
boxplot(DeltaSVM, names=Nb_simul_Sub, notch=T, outline=F, las=1, cex.axis=1.4, cex.lab=1.4,
        xlab="substitution number", ylab="delta SVM")
abline(v=9.5, lty=3)

nbSub <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$NbSub)
boxplot(nbSub, names=Nb_simul_Sub, notch=T, outline=F, las=1,cex.axis=1.4, cex.lab=1.4,
        xlab="Nb mutations steps", ylab="substitution number")

abline(h=c(10, 20, 40, 60, 80, 100), lty=3)
abline(v=9.5, lty=3)

propSignif <- lapply(Nb_simul_Sub, function(x) length(which(Simulations[[x]]$FDR<0.1))*100/nrow(Simulations[[x]]))
boxplot(propSignif, names=Nb_simul_Sub,cex.axis=1.4, cex.lab=1.4,
        xlab="substitution number", ylab="% positive TFBS")
abline(v=9.5, lty=3)
abline(h=c(20, 40, 60, 80), lty=3)

pval <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$pval.high)
boxplot(pval, names=Nb_simul_Sub, notch=T, outline=F, xlab="Nb mutations steps", ylab="pval high")

FDR <- lapply(Nb_simul_Sub, function(x) Simulations[[x]]$FDR)
boxplot(FDR, names=Nb_simul_Sub, notch=F, outline=F, xlab="Nb mutations steps", ylab="FDR")


hist(Simulations[[4]]$deltaSVM, breaks=50, xlab="deltaSVM", main="3 mutational steps", xlim=c(-25, 40))
hist(Simulations[[4]]$med.deltaSVM.simul, breaks=50, xlab="random deltaSVM", main="3 mutational steps", xlim=c(-25, 40))

hist(Simulations[[20]]$deltaSVM, breaks=100, xlab="deltaSVM", main="10 mutational steps", xlim=c(-25, 40))
hist(Simulations[[20]]$med.deltaSVM.simul, breaks=50, xlab="random deltaSVM", main="10 mutational steps", xlim=c(-25, 40))
dev.off()

##############################################################################
# Simulation oriented delta
Nb_simul_Sub <- c(2, 3, 4, 10, 20, 50)

pdf(paste0(path, "/figures/human_simulations_oriented.pdf"), width=8)
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

dev.off()


random_delta3 <- read.table(paste0(path, "positive_selection/human/simulation_mutational_steps/Distrib_3_mutations.txt"),h=F)
random_delta20 <- read.table(paste0(path, "positive_selection/human/simulation_mutational_steps/Distrib_20_mutations.txt"),h=F)
random_delta40 <- read.table(paste0(path, "positive_selection/human/simulation_mutational_steps/Distrib_40_mutations.txt"),h=F)
random_delta80 <- read.table(paste0(path, "positive_selection/human/simulation_mutational_steps/Distrib_80_mutations.txt"),h=F)

par(mfrow=(c(2,2)))
par(mai=c(0.6,0.8,0.4,0.3), mgp=c(2.2,0.8,0))
simul=as.numeric(random_delta3$V1)
hist(simul, breaks=8000, xlim=c(-30, 30), freq=F, col=rgb(0, 0, 255, alpha = 125, max=255),
     main="3 substitutions", xlab="deltaSVM")
hist(Simulations[[3]]$deltaSVM, breaks=50, add=T, col=rgb(255, 0, 0, alpha = 125, max=255), freq=F)
abline(v=median(simul, na.rm = T), col="blue")
abline(v=median(Simulations[[3]]$deltaSVM), col="red")

simul=as.numeric(random_delta20$V1)
hist(simul, breaks=400, xlim=c(-50, 50), freq=F, col=rgb(0, 0, 255, alpha = 125, max=255),
     main="20 substitutions", xlab="deltaSVM")
hist(Simulations[[20]]$deltaSVM, breaks=100, add=T, col=rgb(255, 0, 0, alpha = 125, max=255), freq=F)
abline(v=median(simul, na.rm = T), col="blue")
abline(v=median(Simulations[[20]]$deltaSVM), col="red")

simul=as.numeric(random_delta40$V1)
hist(simul, breaks=10000, xlim=c(-50, 50), freq=F, col=rgb(0, 0, 255, alpha = 125, max=255),
     main="40 substitutions", xlab="deltaSVM")
hist(Simulations[[40]]$deltaSVM, breaks=100, add=T, col=rgb(255, 0, 0, alpha = 125, max=255), freq=F)
abline(v=median(simul, na.rm = T), col="blue")
abline(v=median(Simulations[[40]]$deltaSVM), col="red")

simul=as.numeric(random_delta80$V1)
hist(simul, breaks=1000, xlim=c(-70, 70), freq=F, col=rgb(0, 0, 255, alpha = 125, max=255),
     main="80 substitutions", xlab="deltaSVM")
hist(Simulations[[80]]$deltaSVM, breaks=100, add=T, col=rgb(255, 0, 0, alpha = 125, max=255), freq=F)
abline(v=median(simul, na.rm = T), col="blue")
abline(v=median(Simulations[[80]]$deltaSVM), col="red")


col=c(rgb(0, 0, 255, alpha = 125, max=255), rgb(255, 0, 0, alpha = 125, max=255))
names(col) = c("simul", "obs")

par(mfrow=(c(2,1)))
par(mai=c(0.8,0.8,0.4,0.3))
simul=rnorm(10000, mean = 0, sd = 3.3)
hist(simul, freq=F, breaks=50, col=col["simul"], main="3 substitutions", xlab="deltaSVM", xlim=c(-30, 30), cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
hist(Simulations[[3]]$deltaSVM, breaks=80, add=T, col=col["obs"], freq=F)
legend("topright", legend=c("neutral", "observed"), fill=col, bty="n", cex=1.3)

simul=rnorm(10000, mean = -8, sd = 6)
hist(simul, freq=F, breaks=50, col=col["simul"],main="8 substitutions", xlab="deltaSVM", xlim=c(-30, 30), cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
hist(Simulations[[8]]$deltaSVM, breaks=80, add=T, col=col["obs"], freq=F)

