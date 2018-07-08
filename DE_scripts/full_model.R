library(DESeq2)
library(RColorBrewer)
library(gplots)
library(tximport)
library(lattice)
source('~/Documents/UCDavis/Whitehead/RNAseq_15killifish/scripts/plotPCAWithSampleNames.R')
source('~/Documents/UCDavis/Whitehead/RNAseq_15killifish/scripts/overLapper_original.R')
counts <- read.csv("/Users/johnsolk/Documents/UCDavis/Whitehead/16killifish_counts_RNAseq_filtered_2July2018_test.csv",stringsAsFactors = FALSE)
counts2 <-read.csv("~/Documents/UCDavis/Whitehead/16killifish_counts_RNAseq_filtered_2July2018.csv",stringsAsFactors = FALSE)
dim(counts)
head(counts)
#design<-read.csv("~/Documents/UCDavis/Whitehead/salinity_killifish_design.csv",header=TRUE)
design <- counts[counts$GeneName == 'Empty',]
design
sp<-as.character(unlist(design[1,]))
sp<-sp[-c(1,2)]
ph<-as.character(unlist(design[2,]))
ph<-ph[-c(1,2)]
cl<-as.character(unlist(design[3,]))
cl<-cl[-c(1,2)]
de<-as.character(unlist(design[4,]))
de<-de[-c(1,2)]

x <- counts2[ -c(1,2) ]
x <- x+1
log_x<-log(x)
colnames(log_x)
names<-colnames(log_x)
names

pca = prcomp(t(log_x))

#fac = factor(sapply(names,function(x){strsplit(x,'.quant')[[1]][1]}))
#fac2 = factor(sapply(fac,function(x){strsplit(x,'_')[[1]][1]}))
#fac= factor(c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt"))
fac = factor(cl)
fac
colours = function(vec){
  #cols=cols=palette(brewer.pal(n=7,name="Dark2"))
  cols=rainbow(length(unique(vec)))
  #print(cols)
  #cols = c('#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')
  #cols = palette(brewer.pal(n=16,name="Dark2"))
  print(cols)
  return(cols[as.numeric(as.factor(vec))])}

#pdf("PCA.pdf", width=11, height=8.5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(pca$x[,1:2], 
     col=colours(fac), 
     #pch=19,
     pch = c(16, 2, 9)[as.numeric(as.factor(ph))],
     cex=2,
     xlab="PC1",
     ylab="PC2",
     cex.lab=2,
     cex.axis = 2)
legend(185,-150,legend=c("Clade 1","Clade 2","Clade 3"),col=rainbow(length(unique(fac))),cex=1.5, pch=19)
legend(80,-150,legend=c("Brackish","Freshwater","Marine"),cex=1.5,pch=c(16, 2, 9))
#dev.off()

geneID <- counts2$GeneName
rownames(counts2) <- counts2$GeneName 
counts2 <- counts2[-c(1)]
colnames(counts2)
  
cols<-colnames(x)
ExpDesign <- data.frame(row.names=cols, physiology = ph,clade = cl)
ExpDesign

all(rownames(ExpDesign) == colnames(counts2))
counts_round<- round(counts2,digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~ clade)
#dds <- DESeqDataSetFromTximport(countData = counts2,colData = ExpDesign,design = ~ clade)
dds<-DESeq(dds,betaPrior=FALSE)
matrix(resultsNames(dds))
log_cds<-rlog(dds)
plotDispEsts(dds)
plotPCAWithSampleNames(log_cds,intgroup="clade",ntop=40000)

counts_table = counts(dds, normalized=TRUE )
dim(counts_table)
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
dim(filtered_norm_counts)
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)


res.1<-results(dds,contrast=c("clade","Clade1","Clade2"))
dim(res.1)
res.2<-results(dds,contrast=c("clade","Clade1","Clade3"))

res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
#res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
#res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="Clade 1 vs. Clade 3 (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")

plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="Clade 1 vs. Clade 2 (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")



counts <- filtered_norm_counts
#a<-seq(1, length(rownames(counts)), by=1)
#rownames(counts)<-a
#row.has.na <- apply(counts, 1, function(x){any(is.na(x))})
#dim(counts)
#final.filtered <- counts[!row.has.na,]
#dim(final.filtered)
#rownames(counts)<-as.numeric(rownames(counts))
drops <- c("GeneID")
counts<-counts[ , !(names(counts) %in% drops)]
colnames(counts)
id <-rownames(counts)
d<-as.matrix(counts)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="16 species killifish, osmotic challenge RNAseq", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)