data<-read.table(file="merged_fcounts.txt", header=T)
data<-data[,-1]
keep<-rowSums(data)>0
dataraw<-data[keep,]
colnames(dataraw)<-substring(colnames(dataraw),11,19)

topGenes<-read.table(file="_filtered.txt", header=F)
colnames(topGenes)<-"genes"
#Up<-subset(dataraw, rownames(dataraw) %in% topGenes$genes)
Up<-dataraw[match(topGenes$genes, rownames(dataraw)), ]

#Normalize dataset
dat <- transform(Up, min = pmin(i701_i505, i703_i502))
dfnorm<-dat/dat$min
dfnorm<-dfnorm[,-13]


pdf("normalized_foldchangediff.pdf")
heatmap.2(data.matrix(dfnorm), Colv=TRUE, Rowv=TRUE, col=bluered, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()

colors = c(seq(0,1, length=100),seq(1.1,5, length=100))
my_palette <- colorRampPalette(c("blue", "red"))

pdf("colorkey_2.pdf")
heatmap.2(data.matrix(dfnorm), Colv=TRUE, Rowv=TRUE, col=my_palette, breaks=colors, scale="none", trace="none", cexCol = 0.8, cexRow=0.4)
dev.off()


