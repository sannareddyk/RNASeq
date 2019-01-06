##data normalization
#datanorm<-rbind(datafilt,tot.counts = colSums(datafilt))
#########
keep<-rowSums(data[,2:13])>0
datafilt<-data[keep,]

##cpm
datacpm<-apply(datafilt[,2:13],2,function(x)(x/sum(x))*1000000)
datacpmgl<-cbind(datafilt$geneLength,datacpm)
colnames(datacpmgl)[1]<-"geneLength"

##rpkm
datacpmgl<-as.data.frame(datacpmgl)
datacpmgl$geneLength<-datacpmgl$geneLength/1000
for (column.name in names(datacpmgl)[2:ncol(datacpmgl)]){datacpmgl[column.name] = datacpmgl[column.name]/datacpmgl$geneLength}
datarpkm<-datacpmgl[,2:13]

logrpkm<-t(log(1+datarpkm,2))
logdata<-na.omit(logrpkm)
scaleddata<-scale(logdata)

##clustering
dist.mat<-dist(logdata, method="euclidean")
dist.mat<-dist(scaleddata, method="euclidean")
clusters<-hclust(dist.mat, method="ward")
plot(clusters)

##pca
pca<-prcomp(logdata)
plot(pca,type="l")
biplot(pca)

pca=prcomp(logdata, center=T, scale=T)
plot(pca$x[, 1], pca$x[, 2])
pca.var<-pca$sdev^2
pca.var.per<-round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main = "Scree plot", xlab ="Principal component", ylab="Percent variation")

pca.data <- data.frame(Sample=rownames(pca$x),
X=pca$x[,1],
Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw()
  ggtitle("logrpkm PCA Lucia")
  
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])
top_10_genes
pca$rotation[top_10_genes,1]

summary(pca)



