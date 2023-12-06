######################################################################

args = commandArgs(trailingOnly=TRUE)
sp = args[1]   

path=paste0("/Users/alaverre/Documents/Detecting_positive_selection/data/genome_sequence/", sp, '/') ## define path $seqens
GTF=list.files(path = path, pattern = "*.gtf")

######################################################################

annot=read.table(paste0(path, GTF), h=F, sep="\t", stringsAsFactors=F, quote="\"")

## split by annotations
TSS=annot[which(annot$V3=="start_codon"),]
TSS$V9 = gsub(";", "", TSS$V9)
TSS$V9 = sapply(strsplit(TSS$V9," "), `[`, 2)

gene=annot[which(annot$V3=="gene"),]
CDS=annot[which(annot$V3=="CDS"),]

non.coding = annot[grep("protein_coding", annot$V9, invert=T),] #exclude coding
non.coding =  non.coding[grep("pseudogene", non.coding$V9, invert=T),] #exclude pseudogene
non.coding.exon=non.coding[which(non.coding$V3=="exon"),]

######################################################################

#write.table(TSS, file=paste(pathGTF, prefixAnnot,"_TSS.gtf",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(gene, file=paste(pathGTF, prefixAnnot,"_genes.gtf",sep=""),row.names=F,col.names=F,sep="\t",quote=F)
write.table(CDS, file=paste0(path, "CDS_", GTF),row.names=F,col.names=F,sep="\t",quote=F)
#write.table(non.coding.exon, file=paste(pathGTF, prefixAnnot,"_noncoding_exons.gtf",sep=""),row.names=F,col.names=F,sep="\t",quote=F)

######################################################################