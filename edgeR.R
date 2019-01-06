library(edgeR)
counts.table="path/to/edgeR/merged_fcounts.txt"
countData<-read.table(counts.table, header = TRUE, row.names = 1)
#countData<-countData[,c("SLX.14534.i701_i505.HMW3WBBXX.s_8","SLX.14534.i703_i502.HMW3WBBXX.s_8","SLX.14534.i702_i505.HMW3WBBXX.s_8","SLX.14534.i704_i502.HMW3WBBXX.s_8")]
#countData<-countData[,c("SLX.14534.i701_i502.HMW3WBBXX.s_8", "SLX.14534.i703_i503.HMW3WBBXX.s_8", "SLX.14534.i702_i502.HMW3WBBXX.s_8", "SLX.14534.i704_i503.HMW3WBBXX.s_8" )]

metadata.file="path/to/metadata.txt"
colData<-read.table(metadata.file, header = FALSE, row.names = 1)
colnames(colData)<-"condition"
d <- DGEList(counts=countData,group=factor(colData$condition))

#cpm cal'n and filtering
keep<-rowSums(cpm(d)>1)>=1
d_filt<-d[keep,,keep.lib.sizes=FALSE]

#Normalization, TMM
d_filt<-calcNormFactors(d_filt)
d_filt$samples

#data exploration
plotMDS(d_filt, method="bcv", col=as.numeric(d_filt$samples$group))
#legend("topright", as.character(unique(d_filt$samples$group)), col=1:2, pch=10)

#dispersion and exactTest
d1 <- estimateCommonDisp(d_filt, verbose=T)
d1<-estimateTagwiseDisp(d1)
names(d1)

de<-exactTest(d1)
topTags(de, n=10)
de1 <- decideTestsDGE(de, adjust.method="BH", p.value=0.05)
summary(de1)
