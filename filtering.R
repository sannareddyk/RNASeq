setwd("")
res<-read.table(file="_results_unfiltered.txt", header=TRUE, row.names = 1)
res_padjfilt<-res[which(res$padj<0.05),]
dim(res_padjfilt)
res_padjlogfcfilt<-res_padjfilt[which(res_padjfilt$log2FoldChange>=2 | res_padjfilt$log2FoldChange<=-2 ),]
write.table(res_padjlogfcfilt, file="resultsfiltered_Pos_vs_Neg.txt", sep="\t", quote=F, row.names=T, col.names=NA )


final_filter<-res_padjlogfcfilt[which(res_padjlogfcfilt$padj<1e-10),]
final_filter<-final_filter[which(final_filter$baseMean>=50),]
final_data<-final_filter[order(final_filter$log2FoldChange),]
write.table(final_data, file="filtered_xx_vsxx.txt", sep="\t", row.names=T, col.names=NA, quote=FALSE)


downreg_genes<-final_data[final_data$log2FoldChange<0,]
downreg_genes<-downreg_genes[order(downreg_genes$log2FoldChange),]
upreg_genes<-final_data[final_data$log2FoldChange>0,]
upreg_genes<-upreg_genes[order(-upreg_genes$log2FoldChange),]

setwd( "path/to/_FASTQ/_featureCounts/")
Lucia_rpkm<-read.table(file="_rpkms.txt", header=T)
Lucia_NB<-Lucia_rpkm[,c("SLX.14534.i701_i505.HMW3WBBXX.s_8","SLX.14534.i703_i502.HMW3WBBXX.s_8","SLX.14534.i702_i505.HMW3WBBXX.s_8","SLX.14534.i704_i502.HMW3WBBXX.s_8")]
keep<-rowSums(Lucia_NB)>0
datarpkm<-Lucia_NB[keep,]
logrpkm<-log(1+datarpkm,2)
Up<-
pdf("logRPKM_Up.pdf")
heatmap.2(data.matrix(Up), Colv=NA, Rowv=TRUE, col=bluered, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()

data<-read.table(file="merged_fcounts.txt", header=T)
data<-data[,-1]
keep<-rowSums(data)>0
dataraw<-data[keep,]
data_NB<-dataraw[,c("SLX.14534.i701_i505.HMW3WBBXX.s_8","SLX.14534.i703_i502.HMW3WBBXX.s_8","SLX.14534.i702_i505.HMW3WBBXX.s_8","SLX.14534.i704_i502.HMW3WBBXX.s_8")]
logdata_NB<-log(1+data_NB,2)
Up<-logdata_NB[rownames(upreg_genes),]
colnames(Up)<-substring(colnames(Up),11,19)
pdf("logCounts_Up.pdf")
heatmap.2(data.matrix(Up), Colv=NA, Rowv=NA, col=bluered, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()


Up<-data_NB[rownames(upreg_genes),]
pdf("rawCounts_Up.pdf")
heatmap.2(data.matrix(dfnorm), Colv=NA, Rowv=TRUE, col=bluered, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()

#Normalize dataset
normalize<- function(x) {
  return (x/min(x))
}
dfnorm<-t(apply(Up,1, normalize))
pdf("normalized_foldchangediff.pdf")
heatmap.2(data.matrix(dfnorm), Colv=NA, Rowv=TRUE, col=bluered, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()

colors = c(seq(1,1.10, length=100),seq(1.11,4, length=100))
my_palette <- colorRampPalette(c("blue", "red"))

pdf("colorkey_9.pdf")
heatmap.2(data.matrix(dfnorm), Colv=NA, Rowv=TRUE, col=my_palette, breaks=colors, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()
