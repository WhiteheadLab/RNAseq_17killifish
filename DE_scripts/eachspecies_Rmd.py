import os
import os.path
import subprocess
from subprocess import Popen, PIPE


def make_header(species):
	header='''
---
title: "RNAseq, DE analysis for {} killifish osmotic challenge"
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
'''.format(species)
	return header

def load_packages():
	packages="""
```{{r LoadPackages, results='hide', include=FALSE}}
packages<-function(x){{
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){{
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }}
}}
bioconductors <- function(x){{
    x<- as.character(match.call()[[2]])
    if (!require(x, character.only = TRUE)){{
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkgs=x)
      require(x, character.only = TRUE)
    }}
}}
packages("RColorBrewer")
packages("pheatmap")
packages("vsn")
packages(pheatmap)
packages(cowplot)
# for DESeq
packages(lattice)
packages(RColorBrewer)
packages(dplyr)
packages(tidyr)
packages(gplots)
packages(ggplot2)
bioconductors(DESeq2)
bioconductors(tximport)
bioconductors(biomaRt)
```
""".format()
	return packages

def get_files():
	files="""
```{{r load custom scripts and data, results='hide', include=FALSE}}
if(!file.exists('plotPCAWithSampleNames.R')){{
  download.file('https://gist.githubusercontent.com/ljcohen/d6cf3367efae60d7fa1ea383fd7a1296/raw/f11030ced8c7953be6fada0325b92c20d369e0a7/plotPCAWithSampleNames.R', 'plotPCAWithSampleNames.R')
	}}
source('plotPCAWithSampleNames.R')
if(!file.exists('overLapper_original.R')){{
  download.file("https://gist.githubusercontent.com/ljcohen/b7e5fa93b8b77c33d26e4c44cadb5bb7/raw/3a2f84be2e937ef721cafaf308ff6456168f30d9/overLapper_original.R",'overLapper_original.R')
	}}
source('overLapper_original.R')
# This is just the counts with Experimental Design Info in the last 5 rows
if(!file.exists('../../Ensembl_species_counts_designfactors.csv')){{
  download.file("https://osf.io/7vp38/download",'../../Ensembl_species_counts_designfactors.csv')
	}}
counts_design <- read.csv("../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)
```
""".format()
	return files

def design(species):
	design_info='''
# Load counts and design:

Dimensions of counts table, and design table:
```{{r format counts and ExpDesign,}}
design <- counts_design[counts_design$Ensembl == 'Empty',]
drops <- c("X","Ensembl")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
design <- design[ , startsWith(names(design),"{}")]
counts <- counts[ , startsWith(names(counts),"{}")]
dim(counts)
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
'''.format(species,species)
	return design_info


def filtering_counts(n_samples,expn):
	filtered_counts="""
# Filtering counts

The following shows the dimensions of the dataframe when we filter out genes with low counts.

{} samples must have a count of at least {}:

```{{r filtering5,}}
filter <- rownames(counts[rowSums(counts >= {}) >= {},])
filtered_counts <- counts[filter,]
dim(filtered_counts)
```
""".format(n_samples,expn,n_samples,expn)
	return filtered_counts

def run_DESeq():
	DESeq_string="""
```{{r run DESeq, results='hide', include=FALSE}}
all(rownames(ExpDesign) == colnames(counts))
counts_round<- round(data.matrix(counts),digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~condition)
dds<-DESeq(dds)
```
""".format()
	return DESeq_string

def DESeq_QC():
	qc="""
# DESeq

Model results, dispersion plot, and mean variance plot.

```{{r DESeq QC, }}
resultsNames(dds)
plotDispEsts(dds)
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
```
""".format()
	return qc

def PCA():
	PCA_plot='''

# PCA

```{{r PCA of counts,}}
plotPCA(vsd, intgroup=c("condition"))
plotPCAWithSampleNames(vsd,intgroup=c("condition"))
```
'''.format()
	return PCA_plot

def norm_counts():
	counts="""
```{{r dim counts, results='hide', include=FALSE}}
# get counts
counts_table = counts(dds, normalized=TRUE)
dim(counts_table)
```

After filtering for low expression (where rowSum is greater than or equal to 1):

```{{r filtering again, results='hide', include=FALSE}}
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
dim(filtered_norm_counts)
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
```
""".format()
	return counts

def MA_plot(contrast1,contrast2,species,n):
	MA="""
# MA plot, {} vs. {}
```{{r MA plot {}, }}
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
res<-results(dds,contrast=c("condition","{}","{}"))
res_ordered <-as.data.frame(res[order(res$padj),])
res_filtered <-subset(res_ordered,res_ordered$padj<0.05)
id<-rownames(res_filtered)
res_filtered<-cbind(res_filtered,id)
plot(log2(res$baseMean), res$log2FoldChange, 
     col=ifelse(res$padj < 0.05, "red","gray67"),
     main="{} ({} vs. {}) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
resSig = res_ordered[rownames(res_ordered) %in% protein_id, ]
dim(resSig)
genes<-rownames(resSig)
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=gene_id,pos=2,cex=0.60)
```
""".format(contrast1,contrast2,n,contrast1,contrast2,species,contrast1,contrast2)
	return MA


def biomaRt():
	query="""
```{{r biomaRt,results='hide', include=FALSE}}
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

# link goi Ensembl ID to external_gene_name or description
gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
```
""".format()
	return query

def gene_plot(n,goi,geneID):
	goi_plot="""
## {}
```{{r plot goi {},}}
tcounts <- t(log2((counts(dds[c("{}"), ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-1+1):ncol(.))

C1<-ggplot(tcounts, aes(condition, expression,group=1)) +
  geom_point() + 
  scale_x_discrete(limits=c('0.2_ppt','transfer','15_ppt')) +
  stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",width=0.2) +
  theme_bw() +
  theme(legend.position="none",panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y="Expression (log2 normalized counts)")+
  ggtitle("{}")
plot(C1)
```
""".format(geneID,n,goi,geneID)

	return goi_plot

def merge_counts_stats_annotations(species):
	merged="""
```{{r merge and write, results='hide', include=FALSE}}
# -----------------------------
# stats results
# -----------------------------

res_BW_v_FW <- results(dds, tidy=TRUE, contrast=c("condition","15_ppt","0.2_ppt")) %>% arrange(padj) %>% tbl_df() 
res_TR_v_FW <- results(dds, tidy=TRUE, contrast=c("condition","transfer","0.2_ppt")) %>% arrange(padj) %>% tbl_df() 
res_TR_v_BW <- results(dds, tidy=TRUE, contrast=c("condition","transfer","15_ppt")) %>% arrange(padj) %>% tbl_df() 

# -----------------------------
# counts 
# -----------------------------

cols <- colnames(counts_table)
counts_table <- as.data.frame(counts_table[,cols])
dim(counts_table)

# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------

names(res_BW_v_FW)[names(res_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_BW_v_FW)[names(res_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------

names(res_TR_v_FW)[names(res_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_TR_v_FW)[names(res_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'

# -----------------------------
# column names for stats from TR_v_BW specific contrast
# -----------------------------

names(res_TR_v_BW)[names(res_TR_v_BW) == 'padj'] <- 'padj-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'baseMean'] <- 'baseMean-ALL'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'lfcSE'] <- 'lfcSE-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'stat'] <- 'stat-TR-v-15ppt'
names(res_TR_v_BW)[names(res_TR_v_BW) == 'pvalue'] <- 'pvalue-TR-v-15ppt'

# -----------------------------
# merge counts and stats
# -----------------------------
res_TR_v_FW <- as.data.frame(res_TR_v_FW)
rownames(res_TR_v_FW) <- res_TR_v_FW$row
counts_table_stats <- merge(as.data.frame(res_TR_v_FW),counts_table,by=0)
counts_table_stats <- merge(as.data.frame(res_BW_v_FW),counts_table_stats,by='row')
counts_table_stats <- merge(as.data.frame(res_TR_v_BW),counts_table_stats,by='row')
dim(counts_table_stats)

# -----------------------------
# merge annotations with stats
# -----------------------------
counts_table_ann <- merge(query,counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
counts_table_ann <- counts_table_ann[!duplicated(counts_table_ann$ensembl_peptide_id), ]
rownames(counts_table_ann) <- counts_table_ann$ensembl_peptide_id
dim(counts_table_ann)
counts_table_ann <- counts_table_ann[ , -which(names(counts_table_ann) %in% c("Row.names"))]
# -----------------------------
# write csv files
# -----------------------------
counts_table_ann <- counts_table_ann[order(counts_table_ann[,19]), ]
write.csv(counts_table_ann,"/Users/johnsolk/Documents/UCDavis/Whitehead/counts_stats_byspecies/{}_stats_annotations_counts.csv")
```
""".format(species)
	return merged 

def subset_sig():
	subset = """

# Significant genes

Number of significant genes between conditions 15ppt vs. 0.2ppt, padj <0.05:
```{{r subset, }}
sig <- subset(counts_table_stats, counts_table_stats$`padj-15ppt-v-0.2ppt`<= 0.05)
dim(sig)
sig_id <- sig$row
counts_table <- counts(dds,normalized=TRUE)
counts_sig <- counts_table[rownames(counts_table) %in% sig_id,]
```
""".format()
	return subset

def heatmap(species):
	heatmap_plot = """
# Heatmap

```{{r heatmap,}}
id <- sig_id
d<-as.matrix(counts_sig)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="{}, padj<0.05",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
```
""".format(species)
	return heatmap_plot

def make_Rmd(outfile,species):
	header=make_header(species)
	packages=load_packages()
	files=get_files()
	design_info=design(species)
	n_samples="2"
	expn="0.1"
	filtered_counts=filtering_counts(n_samples,expn)
	DESeq_string=run_DESeq()
	qc=DESeq_QC()
	PCA_plot=PCA()
	MA1=MA_plot("15_ppt","0.2_ppt",species,"1")
	MA2=MA_plot("transfer","0.2_ppt",species,"2")
	MA3=MA_plot("transfer","15_ppt",species,"3")
	chunks1=[header,packages,files,design_info,filtered_counts,DESeq_string,qc,PCA_plot,MA1,MA2,MA3]
	counts=norm_counts()
	query=biomaRt()
	merged=merge_counts_stats_annotations(species)
	subset = subset_sig()
	heatmap_plot = heatmap(species)
	chunks2=[counts,query,merged,subset,heatmap_plot]
	genes_proteins = {"ENSFHEP00000000036":"avpr2aa",
	"ENSFHEP00000001609":"slc24a5",
	"ENSFHEP00000003908":"CLDN4",
	"ENSFHEP00000006725":"aqp3",
	"ENSFHEP00000008393":"cftr",
	"ENSFHEP00000009753":"kcnj2a",
	"ENSFHEP00000013324":"polyamine-modulated factor 1-like",
	"ENSFHEP00000015383":"kcnj1a.6",
	"ENSFHEP00000015765":"sept2B",
	"ENSFHEP00000016853":"septin-2", 
	"ENSFHEP00000017303":"cipcb",
	"ENSFHEP00000019510":"clcn2c",
	"ENSFHEP00000025841":"zymogen granule membrane protein 16",
	"ENSFHEP00000031108":"atp1a1b",
	"ENSFHEP00000034177":"solute carrier family 24 member 2"}
	with open(outfile,"w") as Rmd:
		for i in chunks1:
			print(i)
			Rmd.write(i + "\n")
		n=0
		Rmd.write("# Salinity-responsive genes of interest"+"\n")
		for goi in genes_proteins.keys():
			n+=1
			geneID = genes_proteins[goi]
			goi_plot=gene_plot(n,goi,geneID)
			Rmd.write(goi_plot + "\n")
		for j in chunks2:
			print(j)
			Rmd.write(j + "\n")
	print("File written:",outfile)
	exec_string='Rscript -e "library(rmarkdown); rmarkdown::render(\'' + outfile + '\')"' 
	print(exec_string)
	s = subprocess.Popen(exec_string, shell=True)
	s.wait()

species_list=["A_xenica","F_catanatus","F_chrysotus","F_diaphanus","F_grandis",
"F_heteroclitusMDPL","F_heteroclitusMDPP","F_notatus","F_olivaceous",
"F_parvapinis","F_rathbuni","F_sciadicus","F_similis","L_goodei","L_parva"]

outdir="../DE_scripts/by_species/"

for species in species_list:
	outfile=outdir+species+"_DE.Rmd"
	make_Rmd(outfile,species)




