path = "/Users/alaverre/Documents/Detecting_positive_selection/results/positive_selection/rat/"

### DATAS
all_IDs = paste(rep(paste0("pos", 1:1000), each=3), c("sub1", "sub2", "sub3"), sep = ":")
all_col = c("seq_name", all_IDs)
deltas <- read.table(paste0(path, "Wilson/CEBPA/all_possible_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = all_col)

maxSub=152 # 101 for mouse 152 for rat
obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
obs <- read.table(paste0(path, "Wilson/CEBPA/observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
row.names(obs) <- obs$seq_name


##### ALL
# Distribution of all substitution fitness
mat <- as.matrix(deltas[,2:3001])
mat_obs <- as.matrix(obs[,5:105])
col=c(rgb(1,0,0,0.7), rgb(0,1,0,0.7))
hist(mat, breaks=300, xlim=c(-20, 20), main="All subsitutions", freq=F, col=col[1], xlab="deltaSVM per mutation", ylim=c(0,0.2))
hist(mat_obs, breaks=300, xlim=c(-20, 20), col=col[2], freq=F, add=T)

legend("topright", legend=c("All sub", "Obs sub"), fill=col, bty="n")
median(mat, na.rm=T)
median(mat_obs, na.rm=T)
abline(v=0, col="red")

boxplot(list(mat, mat_obs), outline=F, notch=T, col=col, names=c("All sub", "Obs sub"), ylab="deltaSVM per mutation")

# Distribution of substitution fitness for the first 20 sequences
par(mfrow=c(1,2))
boxplot(t(mat[1:20,]), outline=F, notch=T, ylab="deltaSVM per mutation", xlab="sequences", main="All substitutions")
boxplot(t(mat_obs[1:20,]), outline=F, ylab="deltaSVM per mutation", xlab="sequences", main="Observed subsitutions")

treshold=1
all_length = apply(deltas[,2:3001], 1, function(x) length(which(!is.na(x))))
all_neg = apply(deltas[,2:3001], 1, function(x) length(which(x <= -treshold)))
all_pos = apply(deltas[,2:3001], 1, function(x) length(which(x >= treshold)))
all_null = apply(deltas[,2:3001], 1, function(x) length(which(x <=treshold & x >= -treshold)))

nb_sub = obs$NbSub
obs_neg = apply(obs[,5:105], 1, function(x) length(which(x <= -treshold)))
obs_pos = apply(obs[,5:105], 1, function(x) length(which(x >= treshold)))
obs_null = apply(obs[,5:105], 1, function(x) length(which(x <=treshold & x >= -treshold)))

ratio <- data.frame(neg=all_neg, null=all_null, pos=all_pos, all_sub=all_length,
                    obs_neg=obs_neg, obs_null=obs_null, obs_pos=obs_pos, obs_sub=nb_sub,
                    row.names = deltas$seq_name)

boxplot(ratio[1:3]/ratio$all_sub, notch=T, names=c("Negative", "Null", "Positive"), ylab="Proportion")

ratio$Prop_obs_neg = ratio$obs_neg/ratio$obs_sub
ratio$Prop_all_neg = ratio$neg/ratio$all_sub
ratio$Prop_obs_null = ratio$obs_null/ratio$obs_sub
ratio$Prop_all_null = ratio$null/ratio$all_sub
ratio$Prop_obs_pos = ratio$obs_pos/ratio$obs_sub
ratio$Prop_all_pos = ratio$null/ratio$all_sub

ratio$Enrichment_neg = ratio$Prop_obs_neg/ratio$Prop_all_neg
ratio$Enrichment_null = ratio$Prop_obs_null/ratio$Prop_all_null
ratio$Enrichment_pos = ratio$Prop_obs_pos/ratio$Prop_all_pos

# Hypergeometric test
ratio$HyperPval_neg = apply(ratio, 1, function(x) phyper(x["obs_neg"]-1, x["neg"], x["all_sub"]-x["neg"], x["obs_sub"], lower.tail=FALSE))
ratio$HyperPval_null = apply(ratio, 1, function(x) phyper(x["obs_null"]-1, x["null"], x["all_sub"]-x["null"], x["obs_sub"], lower.tail=FALSE))
ratio$HyperPval_pos = apply(ratio, 1, function(x) phyper(x["obs_pos"]-1, x["pos"], x["all_sub"]-x["pos"], x["obs_sub"], lower.tail=FALSE))

# Binomial test
ratio$BinomPval_neg = apply(ratio, 1, function(x) binom.test(x["obs_neg"], x["obs_sub"], x["neg"]/x["all_sub"], alternative = "greater")$p.value)
ratio$BinomPval_null = apply(ratio, 1, function(x) binom.test(x["obs_null"], x["obs_sub"], x["null"]/x["all_sub"], alternative = "greater")$p.value)
ratio$BinomPval_pos = apply(ratio, 1, function(x) binom.test(x["obs_pos"], x["obs_sub"], x["pos"]/x["all_sub"],  alternative = "greater")$p.value)

#hist(ratio$r_pos, breaks=100)
boxplot(list(ratio$Enrichment_neg, ratio$Enrichment_null, ratio$Enrichment_pos), names=c("Negative", "Null", "Positive"), ylab="Ratio", main=paste("Treshold =", treshold))


## Compare Tests
test <- read.table(paste0(path, "/CEBPA_PosSelTest_deltaSVM_10000permutations.txt"), h=T)
row.names(test) <- test$ID
test <- test[rownames(obs),]

# Delta vs TestPos
plot(obs$deltaSVM, test$pval.high, cex=0.1, xlab="delta SVM", ylab="TestPos pval")
cor.test(obs$deltaSVM, test$pval.high)

# Delta vs Enrich pval
plot(obs$deltaSVM, ratio[rownames(obs),"BinomPval_pos"], cex=0.1, xlab="deltaSVM", ylab="EnrichPos pval")
cor.test(obs$deltaSVM, ratio[rownames(obs),"BinomPval_pos"])

# TestPos vs Enrich pval
plot(test$pval.high, ratio[rownames(obs),"BinomPval_pos"], cex=0.1, xlab="TestPos pval", ylab="EnrichPos pval")
cor.test(test$pval.high, ratio[rownames(obs),"BinomPval_pos"])

# TestPos vs NbSub
plot(test$pval.high, test$NbSub, cex=0.1, xlab="TestPos pval", ylab="Nb substitutions", ylim=c(0,30))
cor.test(test$pval.high, test$NbSub)

# Enrich vs NbSub
plot(ratio[rownames(obs),"BinomPval_pos"], test$NbSub, cex=0.1, xlab="EnrichPos pval", ylab="Nb substitutions", ylim=c(0,30))
cor.test(ratio[rownames(obs),"BinomPval_pos"], test$NbSub)

