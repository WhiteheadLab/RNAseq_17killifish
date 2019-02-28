



This report was automatically generated with the R package **knitr**
(version 1.21).


```r

---
title: "RNAseq, DE analysis for A_xenica killifish osmotic challenge"
author: "Lisa K. Johnson"
date: '`r Sys.Date()`'
output:
  html_document:
    code_folding: hide
    collapsed: no
    df_print: paged
    number_sections: yes
    theme: cerulean
    toc: yes
    toc_depth: 5
    toc_float: yes
---

```

```
## Error: <text>:19:0: unexpected end of input
## 17: 
## 18: 
##    ^
```







```r
design <- counts_design[counts_design$Ensembl == 'Empty',]
drops <- c("X","Ensembl")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
design <- design[ , startsWith(names(design),"A_xenica")]
counts <- counts[ , startsWith(names(counts),"A_xenica")]
dim(design)
```

```
## [1] 5 9
```

```r
dim(counts)
```

```
## [1] 30466     9
```

```r
# design cateogories (full)
species<-as.character(unlist(design[1,]))
nativephysiology<-as.character(unlist(design[2,]))
clade<-as.character(unlist(design[3,]))
condition<-as.character(unlist(design[5,]))
cols<-colnames(counts)
ExpDesign <- data.frame(row.names=cols,
                        condition=condition)
ExpDesign
```

```
##                           condition
## A_xenica_BW_1.quant          15_ppt
## A_xenica_BW_2.quant          15_ppt
## A_xenica_BW_3.quant          15_ppt
## A_xenica_FW_1.quant         0.2_ppt
## A_xenica_FW_2.quant         0.2_ppt
## A_xenica_FW_3.quant         0.2_ppt
## A_xenica_transfer_1.quant  transfer
## A_xenica_transfer_2.quant  transfer
## A_xenica_transfer_3.quant  transfer
```


# Filtering counts

The following shows the dimensions of the dataframe when we filter out genes with low counts.

2 samples must have a count of at least 0.1:


```r
filter <- rownames(counts[rowSums(counts >= 2) >= 0.1,])
filtered_counts <- counts[filter,]
dim(filtered_counts)
```

```
## [1] 21251     9
```





# DESeq


```r
plotDispEsts(dds)
```

<img src="figure/A-xenica-DE-RmdDESeq QC-1.png" title="plot of chunk DESeq QC" alt="plot of chunk DESeq QC" style="display: block; margin: auto;" />

```r
resultsNames(dds)
```

```
## [1] "Intercept"                     "condition_15_ppt_vs_0.2_ppt"  
## [3] "condition_transfer_vs_0.2_ppt"
```

```r
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
```

<img src="figure/A-xenica-DE-RmdDESeq QC-2.png" title="plot of chunk DESeq QC" alt="plot of chunk DESeq QC" style="display: block; margin: auto;" />



# PCA


```r
plotPCA(vsd, intgroup=c("condition"))
```

<img src="figure/A-xenica-DE-RmdPCA of counts-1.png" title="plot of chunk PCA of counts" alt="plot of chunk PCA of counts" style="display: block; margin: auto;" />

```r
plotPCAWithSampleNames(vsd,intgroup=c("condition"))
```

```
## 
## Attaching package: 'genefilter'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     rowSds, rowVars
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4    PC5     PC6     PC7
## Standard deviation     9.7726 7.3272 6.6912 6.1279 5.9461 5.66214 5.24772
## Proportion of Variance 0.2713 0.1525 0.1272 0.1067 0.1004 0.09106 0.07822
## Cumulative Proportion  0.2713 0.4238 0.5510 0.6576 0.7580 0.84910 0.92732
##                            PC8       PC9
## Standard deviation     5.05826 1.262e-14
## Proportion of Variance 0.07268 0.000e+00
## Cumulative Proportion  1.00000 1.000e+00
```

<img src="figure/A-xenica-DE-RmdPCA of counts-2.png" title="plot of chunk PCA of counts" alt="plot of chunk PCA of counts" style="display: block; margin: auto;" />


# MA plot, 15_ppt vs. 0.2_ppt

```r
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
res<-results(dds,contrast=c("condition","15_ppt","0.2_ppt"))
res_ordered <-as.data.frame(res[order(res$padj),])
res_filtered <-subset(res_ordered,res_ordered$padj<0.05)
id<-rownames(res_filtered)
res_filtered<-cbind(res_filtered,id)
plot(log2(res$baseMean), res$log2FoldChange, 
     col=ifelse(res$padj < 0.05, "red","gray67"),
     main="A_xenica (15_ppt vs. 0.2_ppt) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
resSig = res_ordered[rownames(res_ordered) %in% protein_id, ]
dim(resSig)
```

```
## [1] 15  6
```

```r
genes<-rownames(resSig)
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=gene_id,pos=2,cex=0.60)
```

<img src="figure/A-xenica-DE-RmdMA plot 1-1.png" title="plot of chunk MA plot 1" alt="plot of chunk MA plot 1" style="display: block; margin: auto;" />


# MA plot, transfer vs. 0.2_ppt

```r
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
res<-results(dds,contrast=c("condition","transfer","0.2_ppt"))
res_ordered <-as.data.frame(res[order(res$padj),])
res_filtered <-subset(res_ordered,res_ordered$padj<0.05)
id<-rownames(res_filtered)
res_filtered<-cbind(res_filtered,id)
plot(log2(res$baseMean), res$log2FoldChange, 
     col=ifelse(res$padj < 0.05, "red","gray67"),
     main="A_xenica (transfer vs. 0.2_ppt) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
resSig = res_ordered[rownames(res_ordered) %in% protein_id, ]
dim(resSig)
```

```
## [1] 15  6
```

```r
genes<-rownames(resSig)
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=gene_id,pos=2,cex=0.60)
```

<img src="figure/A-xenica-DE-RmdMA plot 2-1.png" title="plot of chunk MA plot 2" alt="plot of chunk MA plot 2" style="display: block; margin: auto;" />


# MA plot, transfer vs. 15_ppt

```r
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
res<-results(dds,contrast=c("condition","transfer","15_ppt"))
res_ordered <-as.data.frame(res[order(res$padj),])
res_filtered <-subset(res_ordered,res_ordered$padj<0.05)
id<-rownames(res_filtered)
res_filtered<-cbind(res_filtered,id)
plot(log2(res$baseMean), res$log2FoldChange, 
     col=ifelse(res$padj < 0.05, "red","gray67"),
     main="A_xenica (transfer vs. 15_ppt) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
resSig = res_ordered[rownames(res_ordered) %in% protein_id, ]
dim(resSig)
```

```
## [1] 15  6
```

```r
genes<-rownames(resSig)
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=gene_id,pos=2,cex=0.60)
```

<img src="figure/A-xenica-DE-RmdMA plot 3-1.png" title="plot of chunk MA plot 3" alt="plot of chunk MA plot 3" style="display: block; margin: auto;" />

# Salinity-responsive genes of interest

## avpr2aa

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000000036"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("avpr2aa")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 1-1.png" title="plot of chunk plot goi 1" alt="plot of chunk plot goi 1" style="display: block; margin: auto;" />


## slc24a5

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000001609"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("slc24a5")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 2-1.png" title="plot of chunk plot goi 2" alt="plot of chunk plot goi 2" style="display: block; margin: auto;" />


## CLDN4

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000003908"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("CLDN4")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 3-1.png" title="plot of chunk plot goi 3" alt="plot of chunk plot goi 3" style="display: block; margin: auto;" />


## aqp3

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000006725"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("aqp3")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 4-1.png" title="plot of chunk plot goi 4" alt="plot of chunk plot goi 4" style="display: block; margin: auto;" />


## cftr

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000008393"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("cftr")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 5-1.png" title="plot of chunk plot goi 5" alt="plot of chunk plot goi 5" style="display: block; margin: auto;" />


## kcnj2a

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000009753"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("kcnj2a")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 6-1.png" title="plot of chunk plot goi 6" alt="plot of chunk plot goi 6" style="display: block; margin: auto;" />


## polyamine-modulated factor 1-like

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000013324"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("polyamine-modulated factor 1-like")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 7-1.png" title="plot of chunk plot goi 7" alt="plot of chunk plot goi 7" style="display: block; margin: auto;" />


## kcnj1a.6

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000015383"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("kcnj1a.6")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 8-1.png" title="plot of chunk plot goi 8" alt="plot of chunk plot goi 8" style="display: block; margin: auto;" />


## sept2B

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000015765"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("sept2B")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 9-1.png" title="plot of chunk plot goi 9" alt="plot of chunk plot goi 9" style="display: block; margin: auto;" />


## septin-2

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000016853"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("septin-2")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 10-1.png" title="plot of chunk plot goi 10" alt="plot of chunk plot goi 10" style="display: block; margin: auto;" />


## cipcb

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000017303"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("cipcb")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 11-1.png" title="plot of chunk plot goi 11" alt="plot of chunk plot goi 11" style="display: block; margin: auto;" />


## clcn2c

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000019510"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("clcn2c")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 12-1.png" title="plot of chunk plot goi 12" alt="plot of chunk plot goi 12" style="display: block; margin: auto;" />


## zymogen granule membrane protein 16

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000025841"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("zymogen granule membrane protein 16")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 13-1.png" title="plot of chunk plot goi 13" alt="plot of chunk plot goi 13" style="display: block; margin: auto;" />


## atp1a1b

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000031108"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("atp1a1b")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 14-1.png" title="plot of chunk plot goi 14" alt="plot of chunk plot goi 14" style="display: block; margin: auto;" />


## solute carrier family 24 member 2

```r
tcounts <- t(log2((counts(dds[c("ENSFHEP00000034177"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression)) +
  geom_point() + 
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("solute carrier family 24 member 2")
plot(C1)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

<img src="figure/A-xenica-DE-Rmdplot goi 15-1.png" title="plot of chunk plot goi 15" alt="plot of chunk plot goi 15" style="display: block; margin: auto;" />



# Get normalized counts and filter out genes with low expression

This is the number of genes in the expression counts table:

```r
# get counts
counts_table = counts(dds, normalized=TRUE)
dim(counts_table)
```

```
## [1] 30466     9
```

After filtering for low expression (where rowSum is greater than or equal to 1):


```r
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
dim(filtered_norm_counts)
```

```
## [1] 16269     9
```

```r
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
```



```r
all_goi<-c("ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
           "ENSFHEP00000015383","ENSFHEP00000009753","ENSFHEP00000006725","ENSFHEP00000008393",
           "ENSFHEP00000013324","ENSFHEP00000001609","ENSFHEP00000013324","ENSFHEP00000034177",
           "ENSFHEP00000015765","ENSFHEP00000017303","ENSFHEP00000000036","ENSFHEP00000031108",
           "ENSFHEP00000016853","ENSFHEP00000003908")

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
ensembl_proteinID = rownames(counts_table)
query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id',
'ensembl_gene_id','gene_biotype','external_gene_name',
'description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
```

```
## Batch submitting query [>-------------------------] 3% eta: 29s
## Batch submitting query [>-------------------------] 5% eta: 26s
## Batch submitting query [=>------------------------] 7% eta: 25s
## Batch submitting query [=>------------------------] 8% eta: 23s
## Batch submitting query [==>-----------------------] 10% eta: 22s
## Batch submitting query [==>-----------------------] 11% eta: 21s
## Batch submitting query [==>-----------------------] 13% eta: 27s
## Batch submitting query [===>----------------------] 15% eta: 26s
## Batch submitting query [===>----------------------] 16% eta: 29s
## Batch submitting query [====>---------------------] 18% eta: 28s
## Batch submitting query [====>---------------------] 20% eta: 29s
## Batch submitting query [=====>--------------------] 21% eta: 28s
## Batch submitting query [=====>--------------------] 23% eta: 26s
## Batch submitting query [=====>--------------------] 25% eta: 25s
## Batch submitting query [======>-------------------] 26% eta: 27s
## Batch submitting query [======>-------------------] 28% eta: 28s
## Batch submitting query [=======>------------------] 30% eta: 27s
## Batch submitting query [=======>------------------] 31% eta: 26s
## Batch submitting query [========>-----------------] 33% eta: 25s
## Batch submitting query [========>-----------------] 34% eta: 23s
## Batch submitting query [========>-----------------] 36% eta: 23s
## Batch submitting query [=========>----------------] 38% eta: 22s
## Batch submitting query [=========>----------------] 39% eta: 21s
## Batch submitting query [==========>---------------] 41% eta: 20s
## Batch submitting query [==========>---------------] 43% eta: 19s
## Batch submitting query [===========>--------------] 44% eta: 18s
## Batch submitting query [===========>--------------] 46% eta: 17s
## Batch submitting query [===========>--------------] 48% eta: 17s
## Batch submitting query [============>-------------] 49% eta: 16s
## Batch submitting query [============>-------------] 51% eta: 15s
## Batch submitting query [=============>------------] 52% eta: 15s
## Batch submitting query [=============>------------] 54% eta: 14s
## Batch submitting query [=============>------------] 56% eta: 13s
## Batch submitting query [==============>-----------] 57% eta: 13s
## Batch submitting query [==============>-----------] 59% eta: 13s
## Batch submitting query [===============>----------] 61% eta: 12s
## Batch submitting query [===============>----------] 62% eta: 12s
## Batch submitting query [================>---------] 64% eta: 11s
## Batch submitting query [================>---------] 66% eta: 10s
## Batch submitting query [================>---------] 67% eta: 10s
## Batch submitting query [=================>--------] 69% eta: 9s
## Batch submitting query [=================>--------] 70% eta: 9s
## Batch submitting query [==================>-------] 72% eta: 8s
## Batch submitting query [==================>-------] 74% eta: 8s
## Batch submitting query [===================>------] 75% eta: 7s
## Batch submitting query [===================>------] 77% eta: 7s
## Batch submitting query [===================>------] 79% eta: 6s
## Batch submitting query [====================>-----] 80% eta: 6s
## Batch submitting query [====================>-----] 82% eta: 5s
## Batch submitting query [=====================>----] 84% eta: 5s
## Batch submitting query [=====================>----] 85% eta: 4s
## Batch submitting query [======================>---] 87% eta: 4s
## Batch submitting query [======================>---] 89% eta: 3s
## Batch submitting query [======================>---] 90% eta: 3s
## Batch submitting query [=======================>--] 92% eta: 2s
## Batch submitting query [=======================>--] 93% eta: 2s
## Batch submitting query [========================>-] 95% eta: 1s Batch
## submitting query [========================>-] 97% eta: 1s Batch submitting
## query [=========================>] 98% eta: 0s Batch submitting query
## [==========================] 100% eta: 0s
```

```r
# link goi Ensembl ID to external_gene_name or description
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
```




```r
# stats results
```

```r
res_BW_v_FW <- results(dds, tidy=TRUE, contrast=c("condition","15_ppt","0.2_ppt")) %>% arrange(padj) %>% tbl_df() 
res_TR_v_FW <- results(dds, tidy=TRUE, contrast=c("condition","transfer","0.2_ppt")) %>% arrange(padj) %>% tbl_df() 
res_TR_v_BW <- results(dds, tidy=TRUE, contrast=c("condition","transfer","15_ppt")) %>% arrange(padj) %>% tbl_df() 
```

```r
# counts 
```

```r
cols <- colnames(counts_table)
counts_table <- as.data.frame(counts_table[,cols])
dim(counts_table)
```

```
## [1] 30466     9
```

```r
# column names for stats from BW_v_FW specific contrast
```

```r
names(res_BW_v_FW)[names(res_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
```

```r
# column names for stats from TR_v_FW specific contrast
```

```r
names(res_TR_v_FW)[names(res_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
```

```r
# column names for stats from TR_v_BW specific contrast
```

```r
names(res_TR_v_BW)[names(res_TR_v_BW) == 'padj'] <- 'padj-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'baseMean'] <- 'baseMean-ALL'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'lfcSE'] <- 'lfcSE-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'stat'] <- 'stat-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'pvalue'] <- 'pvalue-TR-v-15ppt'
```

```r
# merge counts and stats
```

```r
res_TR_v_FW <- as.data.frame(res_TR_v_FW)
rownames(res_TR_v_FW) <- res_TR_v_FW$row
counts_table_stats <- merge(as.data.frame(res_TR_v_FW),counts_table,by=0)
counts_table_stats <- merge(as.data.frame(res_BW_v_FW),counts_table_stats,by='row')
counts_table_stats <- merge(as.data.frame(res_TR_v_BW),counts_table_stats,by='row')
dim(counts_table_stats)
```

```
## [1] 30466    29
```

```r
# merge annotations with stats
```

```r
counts_table_ann <- merge(query,counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
counts_table_ann <- counts_table_ann[!duplicated(counts_table_ann$ensembl_peptide_id), ]
rownames(counts_table_ann) <- counts_table_ann$ensembl_peptide_id
dim(counts_table_ann)
```

```
## [1] 30466    35
```

```r
counts_table_ann <- counts_table_ann[ , -which(names(counts_table_ann) %in% c("Row.names"))]
```

```r
# write csv files
```

```r
counts_table_ann <- counts_table_ann[order(counts_table_ann[,19]), ]
write.csv(counts_table_ann,"/Users/johnsolk/Documents/UCDavis/Whitehead/counts_stats_byspecies/A_xenica_stats_annotations_counts.csv")
```

```

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.3
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] genefilter_1.64.0           hexbin_1.27.2              
##  [3] biomaRt_2.38.0              tximport_1.10.1            
##  [5] DESeq2_1.22.2               SummarizedExperiment_1.12.0
##  [7] DelayedArray_0.8.0          BiocParallel_1.16.6        
##  [9] matrixStats_0.54.0          GenomicRanges_1.34.0       
## [11] GenomeInfoDb_1.18.2         IRanges_2.16.0             
## [13] S4Vectors_0.20.1            gplots_3.0.1.1             
## [15] tidyr_0.8.2                 dplyr_0.8.0.1              
## [17] lattice_0.20-38             cowplot_0.9.4              
## [19] ggplot2_3.1.0               vsn_3.50.0                 
## [21] Biobase_2.42.0              BiocGenerics_0.28.0        
## [23] pheatmap_1.0.12             RColorBrewer_1.1-2         
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6           bit64_0.9-7            httr_1.4.0            
##  [4] progress_1.2.0         tools_3.5.2            backports_1.1.3       
##  [7] R6_2.4.0               affyio_1.52.0          rpart_4.1-13          
## [10] KernSmooth_2.23-15     Hmisc_4.2-0            DBI_1.0.0             
## [13] lazyeval_0.2.1         colorspace_1.4-0       nnet_7.3-12           
## [16] withr_2.1.2            prettyunits_1.0.2      tidyselect_0.2.5      
## [19] gridExtra_2.3          curl_3.3               bit_1.1-14            
## [22] compiler_3.5.2         preprocessCore_1.44.0  htmlTable_1.13.1      
## [25] labeling_0.3           caTools_1.17.1.1       scales_1.0.0          
## [28] checkmate_1.9.1        affy_1.60.0            stringr_1.4.0         
## [31] digest_0.6.18          foreign_0.8-71         XVector_0.22.0        
## [34] base64enc_0.1-3        pkgconfig_2.0.2        htmltools_0.3.6       
## [37] highr_0.7              limma_3.38.3           htmlwidgets_1.3       
## [40] rlang_0.3.1            rstudioapi_0.9.0       RSQLite_2.1.1         
## [43] gtools_3.8.1           acepack_1.4.1          RCurl_1.95-4.11       
## [46] magrittr_1.5           GenomeInfoDbData_1.2.0 Formula_1.2-3         
## [49] Matrix_1.2-15          Rcpp_1.0.0             munsell_0.5.0         
## [52] stringi_1.3.1          zlibbioc_1.28.0        plyr_1.8.4            
## [55] blob_1.1.1             grid_3.5.2             gdata_2.18.0          
## [58] crayon_1.3.4           splines_3.5.2          annotate_1.60.0       
## [61] hms_0.4.2              locfit_1.5-9.1         knitr_1.21            
## [64] pillar_1.3.1           geneplotter_1.60.0     XML_3.98-1.17         
## [67] glue_1.3.0             evaluate_0.13          latticeExtra_0.6-28   
## [70] data.table_1.12.0      BiocManager_1.30.4     gtable_0.2.0          
## [73] purrr_0.3.0            assertthat_0.2.0       xfun_0.5              
## [76] xtable_1.8-3           survival_2.43-3        tibble_2.0.1          
## [79] memoise_1.1.0          AnnotationDbi_1.44.0   cluster_2.0.7-1
```

```r
Sys.time()
```

```
## [1] "2019-02-27 15:27:55 PST"
```

