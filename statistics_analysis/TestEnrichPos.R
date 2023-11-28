path = "/Users/alaverre/Documents/Detecting_positive_selection/results/"
species = c("human", "caroli", "spretus", "mouse", "dog", "rat", "macaca")
TFS = c("CEBPA", "FOXA1", "HNF4A", "HNF6")
treshold=0.5
maxSub=150

####################################################################################################################
output_all_enrichment = paste0(path, "positive_selection/all_deltas/all_Enrichment_treshold_", treshold, ".Rds")

if (file.exists(output_all_enrichment)){
  print("Already done!")
  allEnrichment <- readRDS(output_all_enrichment)
  
  }else{
    
    allEnrichment <- list()
    pdf(paste0(path, "figures/proportion_sub_fitness_class_treshold", treshold, ".pdf"))
    par(mfrow=c(2,2))
    for (sp in species){
      for (TF in TFS){
        
        file_all_deltas = paste0(path, "positive_selection/all_deltas/", sp, "/", TF, "/all_possible_deltaSVM.txt")
        if (file.exists(file_all_deltas)){
          message(sp, " ", TF)
          
          ### DATAS
          all_IDs = paste(rep(paste0("pos", 1:1000), each=3), c("sub1", "sub2", "sub3"), sep = ":") # because max seq length=1000
          all_col = c("seq_name", all_IDs)
          deltas <- read.table(file_all_deltas, h=F, sep="\t", quote="", fill=T, col.names = all_col)
          row.names(deltas) <- deltas$seq_name
          message("Nb sequences: ", nrow(deltas))
          
          obs_col = c("seq_name", "SVM", "deltaSVM", "NbSub", paste("sub", 1:maxSub, sep = ":"))
          obs <- read.table(paste0(path, "positive_selection/all_deltas/", sp, "/", TF, "/observed_deltaSVM.txt"), h=F, sep="\t", quote="", fill=T, col.names = obs_col)
          row.names(obs) <- obs$seq_name
          obs <- obs[rownames(deltas),]
          
          ### Extract only deltas
          deltas_values <- deltas[,2:3001]
          obs_values <- obs[,5:105]
          
          all_length = apply(deltas_values, 1, function(x) length(which(!is.na(x))))
          all_neg = apply(deltas_values, 1, function(x) length(which(x <= -treshold)))
          all_pos = apply(deltas_values, 1, function(x) length(which(x >= treshold)))
          all_null = apply(deltas_values, 1, function(x) length(which(x <=treshold & x >= -treshold)))
          
          nb_sub = obs$NbSub
          obs_neg = apply(obs_values, 1, function(x) length(which(x <= -treshold)))
          obs_pos = apply(obs_values, 1, function(x) length(which(x >= treshold)))
          obs_null = apply(obs_values, 1, function(x) length(which(x <=treshold & x >= -treshold)))
          
          ratio <- data.frame(neg=all_neg, null=all_null, pos=all_pos, all_sub=all_length,
                              obs_neg=obs_neg, obs_null=obs_null, obs_pos=obs_pos, obs_sub=nb_sub,
                              row.names = deltas$seq_name)
          
          boxplot(ratio[1:3]/ratio$all_sub, notch=T, names=c("Negative", "Null", "Positive"), ylab="Proportion", main=paste(sp, TF))
          
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
          #ratio$HyperPval_neg = apply(ratio, 1, function(x) phyper(x["obs_neg"]-1, x["neg"], x["all_sub"]-x["neg"], x["obs_sub"], lower.tail=FALSE))
          #ratio$HyperPval_null = apply(ratio, 1, function(x) phyper(x["obs_null"]-1, x["null"], x["all_sub"]-x["null"], x["obs_sub"], lower.tail=FALSE))
          #ratio$HyperPval_pos = apply(ratio, 1, function(x) phyper(x["obs_pos"]-1, x["pos"], x["all_sub"]-x["pos"], x["obs_sub"], lower.tail=FALSE))
          
          # Binomial test
          ratio$BinomPval_neg = apply(ratio, 1, function(x) binom.test(x["obs_neg"], x["obs_sub"], x["neg"]/x["all_sub"], alternative = "greater")$p.value)
          ratio$BinomPval_null = apply(ratio, 1, function(x) binom.test(x["obs_null"], x["obs_sub"], x["null"]/x["all_sub"], alternative = "greater")$p.value)
          ratio$BinomPval_pos = apply(ratio, 1, function(x) binom.test(x["obs_pos"], x["obs_sub"], x["pos"]/x["all_sub"],  alternative = "greater")$p.value)
          
          message("Prop EnrichPos: ", length(which(ratio$BinomPval_pos<0.05))/nrow(ratio))
          allEnrichment[[sp]][[TF]] <- ratio
        }
      }
    }
    dev.off()
    saveRDS(allEnrichment, file=output_all_enrichment)
  
}

signif <- lapply(allEnrichment, function(x) lapply(x, function(y) length(which(y$BinomPval_pos<0.05))/nrow(y)))
color=c(rep("red", 4), rep("green", 2), rep("yellow", 1), rep("grey", 1), rep("blue", 4), rep("pink", 2), rep("orange", 4))

pdf(paste0(path, "/figures/proportion_EnrichPos_Binom_0.05_treshold_", treshold, ".pdf"))
plot(unlist(signif), col=color, pch=19, xlab="samples", ylab="Proportion EnrichPos")
legend("topleft", legend=species, fill=unique(color), bty="n", ncol=2)
dev.off()

####################################################################################################################
#hist(ratio$r_pos, breaks=100)
boxplot(list(ratio$Enrichment_neg, ratio$Enrichment_null, ratio$Enrichment_pos), names=c("Negative", "Null", "Positive"), ylab="Ratio", main=paste("Treshold =", treshold))

## Compare Tests
pdf(paste0(path, "/figures/Tests_comparison_", sp, "_", TF, "_treshold_", treshold, ".pdf"))
test <- read.table(paste0(path, "positive_selection/", sp, "/", TF, "_PosSelTest_deltaSVM_10000permutations.txt"), h=T)
row.names(test) <- test$ID
test <- test[rownames(obs),]

# Delta vs TestPos
plot(obs$deltaSVM, test$pval.high, cex=0.1, xlab="delta SVM", ylab="Permutations pval", main="Test DeltaSVM vs permutations")
cor <- cor.test(obs$deltaSVM, test$pval.high)
mtext(paste("R2=", signif(cor$estimate,3)))

# Delta vs Enrich pval
plot(obs$deltaSVM, ratio[rownames(obs),"BinomPval_pos"], cex=0.1, xlab="deltaSVM", ylab="EnrichPos pval", main="Test DeltaSVM vs Enrich Pos Sub")
cor <- cor.test(obs$deltaSVM, ratio[rownames(obs),"BinomPval_pos"])
mtext(paste("R2=", signif(cor$estimate,3)))

# TestPos vs Enrich pval
plot(test$pval.high, ratio[rownames(obs),"BinomPval_pos"], cex=0.1, xlab="Permutations pval", ylab="EnrichPos pval", main="Comparison of tests")
cor <- cor.test(test$pval.high, ratio[rownames(obs),"BinomPval_pos"])
mtext(paste("R2=", signif(cor$estimate,3)))

# TestPos vs NbSub
test$Sub_class <- cut(test$NbSub, breaks=c(2, 3, 5, 10, max(test$NbSub)))
boxplot(test$pval.high~test$Sub_class, cex=0.1, ylab="Permutations pval", xlab="Nb substitutions", main="Permutations Test", notch=T)
cor <- cor.test(test$pval.high, test$NbSub)
mtext(paste("R2=", signif(cor$estimate,3)))

# Enrich vs NbSub
boxplot(ratio[rownames(obs),"BinomPval_pos"]~test$Sub_class, cex=0.1, ylab="EnrichPos pval", xlab="Nb substitutions", main="EnrichPos Test", notch=T)
cor <- cor.test(ratio[rownames(obs),"BinomPval_pos"], test$NbSub)
mtext(paste("R2=", signif(cor$estimate,3)))

dev.off()
####################################################################################################################