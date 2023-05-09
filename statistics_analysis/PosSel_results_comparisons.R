##########################################################################################
species = "human" #"mouse" or "human"
sample = "CEBPA" #"CEBPA" #"HNF4A"

path = "/Users/alaverre/Documents/Detecting_positive_selection/"
#pdf(file=paste0(path, "Comparisons_gkm.pdf"))

prefix = ifelse(species == "human", "hsap_", "bl6_")

pathData = paste0(path, "/Tools/JialinTool/data/", species, "/")
pathResults = paste0(path, "/results/tools_tests/", species, "/", sample, "/Test-models/")

##########################################################################################
###### Compare predicted weights between runs ###### 
# Jialin score
predicted.score = read.table(paste0(pathData, species, "_SVM_model/", sample, "/kmer_10_library_weigths.txt"),
                          row.names = 1, col.names = c("kmer", "Jialin"))

# Predicted score from several runs 
cor <- c()
cor_run1 <- c()
for (i in seq(1,10)){
  run = read.table(paste0(pathResults, "/run_", i, "/prediction.txt"), row.names = 1)
  predicted.score[[paste0("run_", i)]] = run$V2
  test = cor.test(predicted.score$Jialin, predicted.score[[paste0("run_", i)]], method="pearson")
  test_run1 = cor.test(predicted.score[["run_1"]], predicted.score[[paste0("run_", i)]], method="pearson")
  cor <- c(cor, test$estimate)
  cor_run1 <- c(cor_run1, test_run1$estimate)
}

# Add run on filtered seauences
run = read.table(paste0(pathResults, "/run_1/filtered.prediction.txt"), row.names = 1)
predicted.score[[paste0("filtered")]] = run$V2
test = cor.test(predicted.score$Jialin, predicted.score[[paste0("filtered")]], method="pearson")
cor <- c(cor, test$estimate)

#test_run1 = cor.test(predicted.score[["run_1"]], predicted.score[[paste0("filtered")]], method="pearson")
#cor_run1 <- c(cor_run1, test_run1$estimate)

#par(mfrow=c(1,2))
boxplot(predicted.score, outline=F, notch=T, main=sample, ylab="Weights")
plot(cor, ylab="Pearson R2 (Jialin vs Run)", xlab="Runs", ylim=c(0.75, 1))

# Add run 1 vs Jialin
test_run1 = cor.test(predicted.score[["run_1"]], predicted.score$Jialin, method="pearson")
cor_run1 <- c(cor_run1, test_run1$estimate)

names(cor_run1) <- c(paste("run_", 1:10), "Jialin")
plot(cor_run1, ylab="Pearson R2 (Run1 vs others)", xlab="Runs", ylim=c(0.75, 1), main=sample)

par(mfrow=c(1,1))
col=c(rgb(0.1,0.4,0.7,0.5),rgb(0.1,0.7,0.4,0.5), rgb(0.7,0.4,0.1,0.5))
boxplot(predicted.score$Jialin, predicted.score$run_1, notch=T, outline=F, names=c("Jialin", "Mine"), col=col, ylab="weight", main=paste(species, sample))
plot(predicted.score$Jialin, predicted.score$run_1, xlab="Jialin weights", ylab="My weights", cex=0.5)
abline(a=0, b=1, col="red")
a=cor.test(predicted.score$Jialin, predicted.score$V2, method="pearson")
mtext(paste0("R2 = ", signif(a$estimate,3)), side=3) 

# Compare DeltaSVM with same ancestral sequences
delta = read.table(paste0(pathData, "/", species, "_deltaSVM/hsap_", sample, "_deltaSVM_highertailTest.txt"))
colnames(delta) = c("ID", "Jialin", "Jialin.SNP", "pval")
rownames(delta) <- delta$ID

# Make comparable ID
splString<-strsplit(delta$ID,"_",fixed=TRUE)
splString<-data.frame(unlist(splString))
ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
ID.bed<-ID.bed[,c(1:3)]
Jialin.ID <- paste(ID.bed[,1], as.numeric(ID.bed[,2])-1, ID.bed[,3], sep=":")

predicted.delta <- delta[,1:2]
SNPnumber <- delta[,c(1,3)]
  
cor <- c()
cor_run1 <- c()
cor_SNP <- c()
for (i in c(seq(1,10))){
  run = read.table(paste0(pathResults, "/run_", i, "/deltaSVM.txt"), row.names = 1, h=T)
  predicted.delta[[paste0("run_", i)]] = run$deltaSVM
  SNPnumber[[paste0("run_", i)]] = run$SNP
  
  test = cor.test(predicted.delta$Jialin, predicted.delta[[paste0("run_", i)]], method="pearson")
  test_run1 = cor.test(predicted.delta[["run_1"]], predicted.delta[[paste0("run_", i)]], method="pearson")
  cor <- c(cor, test$estimate)
  cor_run1 <- c(cor_run1, test_run1$estimate)
  
  test_SNP = cor.test(SNPnumber$Jialin.SNP, SNPnumber[[paste0("run_", i)]], method="pearson")
  cor_SNP <- c(cor_SNP, test_SNP$estimate)
}

# Add run 1 vs Jialin
test_run1 = cor.test(predicted.delta[["run_1"]], predicted.delta$Jialin, method="pearson")
cor_run1 <- c(cor_run1, test_run1$estimate)

predicted.delta <- predicted.delta[,2:ncol(predicted.delta)]
SNPnumber <- SNPnumber[,2:ncol(predicted.delta)]

# Compare DeltaSVM with complete new pipeline
new.pipeline = read.table(paste0(path, "/results/human/CEBPA/PosSelTest_deltaSVM.1000.txt"), h=T, row.names = 1)
new.pipeline = new.pipeline[Jialin.ID,]
predicted.delta[["new.pipe"]] = new.pipeline$deltaSVM
SNPnumber[["new.pipe"]] = new.pipeline$SNP
test=cor.test(predicted.delta$Jialin, new.pipeline$deltaSVM, method="pearson")
test_run1 = cor.test(predicted.delta[["run_1"]], new.pipeline$deltaSVM, method="pearson")
test_SNP = cor.test(SNPnumber$Jialin.SNP, new.pipeline$SNP, method="pearson")
cor <- c(cor, test$estimate)
cor_run1 <- c(cor_run1, test_run1$estimate)
cor_SNP <- c(cor_SNP, test_SNP$estimate)

names(cor_run1) <- c(paste0("run_", 1:10), "Jialin", "New")

boxplot(predicted.delta, outline=F, notch=T, ylab="Delta SVM", main=sample)
plot(cor, ylab="Pearson R2 delta SVM (Jialin vs Run)", xlab="Runs", ylim=c(0.75, 1), main=sample)
plot(cor_run1, ylab="Pearson R2 delta SVM (Run1 vs others)", xlab="Runs", ylim=c(0.9, 1),  main=sample)

par(mfrow=c(1,3))
plot(predicted.delta$run_1~predicted.delta$Jialin, xlab="Jialin deltaSVM", ylab="Run 1 deltaSVM", cex=0.5, cex.lab=1.3, main=paste0("R2=", signif(cor_run1["Jialin"])))
abline(a=0, b=1, col="red")

plot(predicted.delta$run_1~predicted.delta$run_2, xlab="Run 2 deltaSVM", ylab="Run 1 deltaSVM", cex=0.5,cex.lab=1.3, main=paste0("R2=", signif(cor_run1["run_2"])))
abline(a=0, b=1, col="red")

plot(predicted.delta$new.pipe~predicted.delta$Jialin, xlab="Jialin deltaSVM", ylab="New deltaSVM", cex=0.5, cex.lab=1.3, main=paste0("R2=", signif(cor_run1["New"])))
abline(a=0, b=1, col="red")

#### Compare number of SNP
par(mfrow=c(1,2))
boxplot(SNPnumber, outline=F, notch=T)
plot(cor_SNP, ylab="Pearson R2 (Jialin vs Run)", xlab="Runs", ylim=c(0.75, 1))

hist(SNPnumber$Jialin.SNP-new.pipeline$SNP, xlab="SNP number difference (Jialin SNP - New inferred SNP)", main=sample, cex.lab=0.9)

delta$deltaSNP <- delta$Jialin.SNP-new.pipeline$SNP
delta$diffdeltaSVM <- delta$Jialin-new.pipeline$deltaSVM

plot(delta$deltaSNP~delta$diffdeltaSVM)

# significant deltaSVM
Jialin.delta.signif = Jialin.delta[which(Jialin.delta$pval<0.05),]
predicted.delta.signif = predicted.delta[which(predicted.delta$pval.high<0.05),]
boxplot(Jialin.delta$deltaSVM, predicted.delta$deltaSVM, notch=T, outline=F, names=c("Jialin", "Mine"), col=col, ylab="signif deltaSVM")



##########################################################################################
### Comparisons of sequences length
my.length <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length.txt"))
my.new.length <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length_new.txt"))
my.corrected.length <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length_corrected.txt"))
Jialin.length <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length_Jialin.txt"))
posSet.length <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length_posSet.txt"))
posmaf.lenght <- read.table(paste0(path, "/results/human/CEBPA/tmp.test/seq_length_pos-maf.txt"))

splString<-strsplit(Jialin.length$V1,"_",fixed=TRUE)
splString<-data.frame(unlist(splString))
ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
ID.bed<-ID.bed[,c(1:3)]

row.names(posSet.length) <- paste(ID.bed[,1], as.numeric(ID.bed[,2])-1, ID.bed[,3], sep=":")
row.names(Jialin.length) <- paste(ID.bed[,1], as.numeric(ID.bed[,2])-1, ID.bed[,3], sep=":")

row.names(my.length) <- my.length$V1
row.names(my.new.length) <- my.new.length$V1
row.names(my.corrected.length) <- my.corrected.length$V1
row.names(posmaf.lenght) <- posmaf.lenght$V1

Jialin.length$posSet <- posSet.length[row.names(Jialin.length),"V2"]
Jialin.length$Mine <- my.length[row.names(Jialin.length),"V2"]
Jialin.length$Mine.new <- my.new.length[row.names(Jialin.length),"V2"]
Jialin.length$Mine.corrected <- my.corrected.length[row.names(Jialin.length),"V2"]
Jialin.length$PosMAF <- posmaf.lenght[row.names(Jialin.length),"V2"]
Jialin.length$BEDlength <- as.numeric(ID.bed[,3])-(as.numeric(ID.bed[,2])-1)

Jialin.length$deltaposSet <- Jialin.length$BEDlength-Jialin.length$posSet
Jialin.length$deltaJialin <- Jialin.length$BEDlength-Jialin.length$V2
Jialin.length$deltaMine <- Jialin.length$BEDlength-Jialin.length$Mine
Jialin.length$deltaMine.new <- Jialin.length$BEDlength-Jialin.length$Mine.new
Jialin.length$deltaMine.corrected <- Jialin.length$BEDlength-Jialin.length$Mine.corrected
Jialin.length$deltaPosMAF <- Jialin.length$BEDlength-Jialin.length$PosMAF

