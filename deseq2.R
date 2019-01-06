setwd("path/to/deseq2")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")

counts.table="path/to/deseq2/merged_fcounts.txt"
metadata.file="path/to/deseq2/metadata.txt"

countData<-read.table(counts.table, header = TRUE, row.names = 1)
colData<-read.table(metadata.file, header = FALSE, row.names = 1)
colnames(colData)<-"condition"
countData<-countData[,match(rownames(colData),colnames(countData))]

dds<-DESeqDataSetFromMatrix(countData, colData = colData, design = ~condition)

#pre-filtering removing genes/rows with 0 or 1 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
#relevel
dds$condition <- relevel(dds$condition, ref="A")

#differential expression analysis using 'DESeq' function
dds <- DESeq(dds, betaPrior=FALSE)
#results function, contrast
res<-results(dds)
#res<-results(dds, independentFiltering=FALSE)

#plot fold change
#plotMA(res, main="DESeq2, ", ylim=c(-4,4))

#export results
resOrdered <- res[order(res$padj),]
#resSig <- subset(resOrdered, padj < 0.1)
write.table(resOrdered, file=".txt", sep="\t", quote=F, row.names=T, col.names=NA)
