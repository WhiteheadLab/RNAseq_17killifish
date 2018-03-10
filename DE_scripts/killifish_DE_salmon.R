library(DESeq2)
library(RColorBrewer)
library(gplots)
library(tximport)
#setwd("~/Documents/UCDavis/Whitehead/osmotic/DE_analysis_by_species/R_scripts_DE_analysis/")
source('~/Documents/UCDavis/Whitehead/plotPCAWithSampleNames.R')
source('~/Documents/UCDavis/Whitehead/overLapper_original.R')
dir="~/Documents/UCDavis/Whitehead/salmon_denovo_by_species"
files_list = list.files("~/Documents/UCDavis/Whitehead/salmon_denovo_by_species/A_xenica")
files <- file.path(dir,"A_xenica",files_list, "quant.sf")
files
print(file.exists(files))
gene_names_annotated_only <- read.csv("~/Documents/UCDavis/Whitehead/annotation_gene_names/A_xenica_gene_names.csv")
gene_names_unique_and_annotated <- read.csv("~/Documents/UCDavis/Whitehead/annotation_gene_names_novel_plus_annotated/A_xenica_gene_names.csv")
colnames(gene_names_unique_and_annotated)
colnames(gene_names_annotated_only)
dim(gene_names_all)
gene_names_annotated_only <- gene_names_annotated_only[,c(2,6)]
gene_names_unique_and_annotated <- gene_names_unique_and_annotated[,c(2,4)]
cols<-c("transcript_id","gene_id")
colnames(gene_names)<-cols
tx2gene<-gene_names_unique_and_annotated
head(tx2gene)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
cols<-files_list
colnames(txi.salmon$counts)<-cols
head(txi.salmon$counts)
dim(txi.salmon$counts)
# Experimental Design
colnames(txi.salmon$counts)
conditions = factor(c("BW","BW","BW","FW","FW","FW","TR","TR","TR"))
ExpDesign <- data.frame(row.names=colnames(txi.salmon$counts), condition = conditions)
ExpDesign
dds <- DESeqDataSetFromTximport(txi.salmon, ExpDesign, ~condition)
dds<-DESeq(dds,betaPrior=FALSE)
matrix(resultsNames(dds))
log_cds<-rlog(dds)
plotDispEsts(dds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)

res.1<-results(dds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(dds,contrast=c("condition","TR","FW"))
res.3<-results(dds,contrast=c("condition","TR","BW"))
resultsNames(dds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
filtered_norm_counts<-res.1[!rowSums(counts_table==0)>=1, ]
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_heteroclitus.MDPL_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)
