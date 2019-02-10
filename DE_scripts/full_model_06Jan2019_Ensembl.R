# for QA
library("RColorBrewer")
library("pheatmap")
library("vsn")
# for DESeq
library(DESeq2)
library(RColorBrewer)
library(tximport)
library(lattice)
library(dplyr)
library(ggplot2)
#library('gtable')
#library('grid')
#library('magrittr')
library(tidyr)
library("biomaRt")
source('~/Documents/UCDavis/Whitehead/RNAseq_15killifish/scripts/plotPCAWithSampleNames.R')
source('~/Documents/UCDavis/Whitehead/RNAseq_15killifish/scripts/overLapper_original.R')
# This is the one with just counts
counts <- read.csv("/Users/johnsolk/Documents/UCDavis/Whitehead/kfish_salmon/Ensembl_species_counts.csv",stringsAsFactors = FALSE)
# This is just the counts with Experimental Design Info in the last 4 rows
counts2 <-read.csv("/Users/johnsolk/Documents/UCDavis/Whitehead/kfish_salmon/Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)
rownames(counts)<-counts$Ensembl
colnames(counts)
drops <- c("Ensembl")
counts <- counts[ , !(names(counts) %in% drops)]
dim(counts)
dim(counts2)
head(counts)
tail(counts)
tail(counts2)
#design<-read.csv("~/Documents/UCDavis/Whitehead/kfish_salmon/OG_combined_replicates",header=FALSE,sep="\t")
design <- counts2[counts2$Ensembl == 'Empty',]
head(design)
dim(design)
design$type <- c("species","native_salinity","clade","group","condition")
colnames(design)
drops <- c("X","Ensembl")
design <- design[ , !(names(design) %in% drops)]
head(design)
dim(design)

#remove transfer treatment
#transfer_samples<-design[design$type=="condition",]
#transfer_samples<-transfer_samples[, transfer_samples[1, ] == c("transfer")]
#transfer_samples<-colnames(transfer_samples)
#BW_FW_counts<-counts[, -which(colnames(counts) %in% transfer_samples)]
#dim(BW_FW_counts)
#rows<-counts$NCBIproteinID
#rownames(BW_FW_counts)<-rows
#proteinID <- BW_FW_counts$NCBIproteinID
#BW_FW_counts<-BW_FW_counts[ -c(1,2,3,4,5) ]
#dim(BW_FW_counts)

# design cateogories (full)
sp<-as.character(unlist(design[1,]))
#sp<-sp[-c(1,2,3,4,5)]
#sp<-sp[-129]
ph<-as.character(unlist(design[2,]))
#ph<-ph[-c(1,2,3,4,5)]
#ph<-ph[-129]
cl<-as.character(unlist(design[3,]))
#cl<-cl[-c(1,2,3,4,5)]
#cl<-cl[-129]
de<-as.character(unlist(design[4,]))
#de<-de[-c(1,2,3,4,5)]
#de<-de[-129]
condition<-as.character(unlist(design[5,]))
#condition<-condition[-c(1,2,3,4,5)]
#condition<-condition[-129]
species_group<-as.vector(paste(de, sp, sep="_"))
species_group
spceis_group<-species_group[-c(129)]
species_condition<-as.vector(paste(sp,condition,sep="_"))
species_condition
species_condition<-species_condition[-c(129)]
species_condition
# normal full counts
x <- counts
# remove transfer samples
#x <-BW_FW_counts

# design categories, remove "transfer" samples
#design_BW_FW <- design[, -which(colnames(design) %in% transfer_samples)]
#dim(design_BW_FW)
#sp<-as.character(unlist(design_BW_FW[1,]))
#sp<-sp[-c(1,2,3,4,5)]
#sp<-sp[-91]
#length(sp)
#ph<-as.character(unlist(design_BW_FW[2,]))
#ph<-ph[-c(1,2,3,4,5)]
#ph<-ph[-91]
#length(ph)
#cl<-as.character(unlist(design_BW_FW[3,]))
#cl<-cl[-c(1,2,3,4,5)]
#cl<-cl[-91]
#length(cl)
#de<-as.character(unlist(design_BW_FW[4,]))
#de<-de[-c(1,2,3,4,5)]
#de<-de[-91]
#length(de)
#condition<-as.character(unlist(design_BW_FW[5,]))
#condition<-condition[-c(1,2,3,4,5)]
#condition<-condition[-91]
#length(condition)
#species_group<-as.vector(paste(de, sp, sep="_"))

# ================================================
# PCA analysis
# ================================================
# Contrast1 = Mean_species_transfer - Mean_species_FW
# Contrast2 = Mean_species_BW - Mean_species_FW
# Contrast3 = Mean_species_transfer - Mean_species_BW

contrast_TR_v_FW <- read.csv("~/Documents/UCDavis/Whitehead/kfish_salmon/TR_BW_contrasts_design.csv")
contrast_BW_v_FW <- read.csv("~/Documents/UCDavis/Whitehead/kfish_salmon/BW_FW_contrasts_design.csv")
contrast_TR_v_BW <- read.csv("~/Documents/UCDavis/Whitehead/kfish_salmon/TR_BW_contrasts_design.csv")

# filter to only those expressed in all genes
# no 0, going across all rows

contrast_TR_v_FW <- contrast_TR_v_FW[-c(29248:29251),]
names<-contrast_TR_v_FW$X
rownames(contrast_TR_v_FW) <- names
head(contrast_TR_v_FW)
dim(contrast_TR_v_FW)
drops <- c("X")
contrast_TR_v_FW <- contrast_TR_v_FW[ , !(names(contrast_TR_v_FW) %in% drops)]
contrast_TR_v_FW <- na.omit(contrast_TR_v_FW)
dim(contrast_TR_v_FW)
df_contrast_TR_v_FW<-data.matrix(contrast_TR_v_FW)
# get rownames where counts are low
filtered_names<-rownames(counts[!rowSums(counts==0)>=1, ])
length(filtered_names)
filtered_contrast_TR_v_FW<- contrast_TR_v_FW[rownames(contrast_TR_v_FW) %in% filtered_names,] 
dim(filtered_contrast_TR_v_FW)


TR_v_FW_design <- contrast_TR_v_FW[c(29248:29251),]  
BW_v_FW_design <- contrast_BW_v_FW[contrast_BW_v_FW$X == '0']
TR_v_BW_design <- contrast_TR_v_BW[c(29248:29251),]

de<-as.character(unlist(TR_v_FW_design[4,]))
de<-de[-c(1)]
ph<-as.character(unlist(TR_v_FW_design[2,]))
ph<-ph[-c(1)]
dim(TR_v_BW_design)
contrast_TR_v_FW <- contrast_TR_v_FW[-c(29248:29251),]
rownames(contrast_TR_v_FW) <- contrast_TR_v_FW$X
drops <- c("X")
contrast_TR_v_FW <- contrast_TR_v_FW[ , !(names(contrast_TR_v_FW) %in% drops)]

x <- filtered_contrast_TR_v_FW
dim(x)
x <- na.omit(x)
dim(x)
x <- data.matrix(x)

#x <- x+1
#log_x<-log(x)
#colnames(log_x)
colnames(x)
#names<-colnames(log_x)
names<-colnames(x)
pca = prcomp(t(x))
#pca = prcomp(t(log_x))



#fac = factor(sapply(names,function(x){strsplit(x,'.quant')[[1]][1]}))
#fac2 = factor(sapply(fac,function(x){strsplit(x,'_')[[1]][1]}))
#fac= factor(c("sal25ppt","sal25ppt","sal25ppt","sal30ppt","sal30ppt","sal30ppt","sal35ppt","sal35ppt"))
fac = factor(ph)
fac
colours = function(vec){
  #cols=palette(brewer.pal(n=3,name="Dark2"))
  cols=rainbow(length(unique(vec)))
  #print(cols)
  #cols = c('#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')
  #cols = palette(brewer.pal(n=16,name="Dark2"))
  print(cols)
  return(cols[as.numeric(as.factor(vec))])}

#pdf("PCA.pdf", width=11, height=8.5)
summary(pca)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(pca$x[,1:2], 
     col=colours(fac), 
     pch=19,
     #pch = c(16, 2, 9)[as.numeric(as.factor(ph))],
     cex=2,
     xlab="PC1",
     ylab="PC2",
     cex.lab=1.5,
     cex.axis = 1)
#cols=palette(brewer.pal(n=3,name="Dark2"))
legend(300000,400000,legend=c("BW","FW","M"),col=rainbow(length(unique(as.factor(ph)))),cex=1.5, pch=19)
#legend(120,-115,legend=c("0.2 ppt","15 ppt"),cex=1.5,pch=c(16, 2, 9))
#text(pca$x[,1:2], labels=names, pos=3)


summary(pca)
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
legend(120,120,legend=c("Clade 1","Clade 2","Clade 3"),col=rainbow(length(unique(fac))),cex=1.5, pch=19)
legend(102,-110,legend=c("Brackish","Freshwater","Marine"),cex=1.5,pch=c(16, 2, 9))
#text(pca$x[,1:2], labels=names, pos=3)
#dev.off()


#==========================================

# DESeq2 analysis

#==========================================

#rownames(counts) <- counts$GeneName 
#counts <- counts[-c(1)]
#counts <- counts[-c(1)]
#colnames(counts)

#cols<-colnames(BW_FW_counts) 
cols<-colnames(counts)
cols
species_group<-species_group[-c(129)]
species_group
condition
condition<-condition[-c(129)]
condition
species<-sp[-c(129)]
species
physiology<-ph[-c(129)]
physiology
clade<-cl[-129]
clade
ExpDesign <- data.frame(row.names=cols, group = species_group,condition=condition, species = species, physiology = physiology,clade=clade,species_condition=species_condition)
ExpDesign
#m1 <- model.matrix(~ species_condition, ExpDesign)
#m2 <- model.matrix(~ species + species_condition, ExpDesign)
m1 <- model.matrix(~ clade + species + condition,ExpDesign)
colnames(m1)
#colnames(m2)
m1
all.zero <- apply(m1, 2, function(x) all(x==0))
all.zero
idx <- which(all.zero)
m1 <- m1[,-idx]
unname(m1)

all(rownames(ExpDesign) == colnames(counts))
#all(rownames(ExpDesign) == colnames(BW_FW_counts))
#counts_round<- round(BW_FW_counts,digits=0)
counts_round<- round(counts,digits=0)
#dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~condition)
#dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~ condition:species)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = m1)
#dds <- DESeqDataSetFromTximport(countData = counts2,colData = ExpDesign,design = ~ clade + physiology + clade:condition)
#dds_speciescondition <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~species_condition)
#dds<-DESeq(dds,betaPrior=FALSE)
dds <- DESeq(dds, full=m1, betaPrior=FALSE)
#dds_speciescondition$condition <- relevel(dds$condition, ref = "untreated")
#dds_speciescondition <- DESeq(dds_speciescondition)
#dds_LRT<-DESeq(dds, test = "LRT", reduced = ~ 1)
matrix(resultsNames(dds))
matrix(resultsNames(dds_speciescondition))


## clean
# if message after the fitting model and testing step: "rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"
# from this reference:
# https://support.bioconductor.org/p/65091/
# filter out low expressing genes
#dds <- estimateSizeFactors(dds)
#nc <- counts(dds, normalized=TRUE)
#filter <- rowSums(nc >= 10) >= 2
#dds <- dds[filter,]
#dds<-DESeq(dds,betaPrior=FALSE)
#dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = ~species_condition)
#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds, maxit=500)
# best to do this:
ddsClean <- dds[which(mcols(dds)$betaConv),]
dds<-ddsClean

#==========================================

# QA of DESeq

#==========================================

 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
#ntd <- normTransform(dds)

#meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
#log_cds<-rlog(dds)
plotDispEsts(dds)


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#log_dds<-rlog(dds)
plotPCA(vsd, intgroup=c("species"))
plotPCAWithSampleNames(vsd,intgroup=c("species"),ntop=40000)

resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-20,20), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)


# ---------------------------
# make a matrix of counts
# ---------------------------


counts_table = counts(dds, normalized=TRUE )
write.csv(counts_table,"~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_results/Ensembl_15killifish_normcounts_nostats.csv")

# ---------------------------

# Get Results

# ---------------------------
resultsNames(dds)
# speciesA_xen is reference level, and condition0.2 ppt is reference level
# Is there an effect of TR condition in A_xen, 0.2ppt?
condition15_ppt
A_xen_results_BW_v_FW <- results(dds, contrast=list(c("condition15_ppt")))
A_xen_results_transfer_v_FW <- results(dds, contrast=list(c("conditiontransfer")))

res <- results(dds, contrast=c("condition","B","A")))


res_FWvBW <- results(dds, tidy=TRUE, contrast=c("species_condition","A_xenica_0.2_ppt","A_xenica_15_ppt")) %>% arrange(padj) %>% tbl_df()
#res_BWvTR <- results(dds, tidy=TRUE, contrast=c("condition", "15_ppt", "transfer")) %>% arrange(padj) %>% tbl_df()
res_FWvTR <- results(dds, tidy=TRUE, contrast=c("condition", "0.2_ppt", "transfer")) %>% arrange(padj) %>% tbl_df()

#res <- results(dds, tidy=TRUE,contrast=list(c("species_condition_A_xenica_15_ppt_vs_A_xenica_0.2_ppt","species_condition_A_xenica_transfer_vs_A_xenica_0.2_ppt","species_condition_F_catanatus_0.2_ppt_vs_A_xenica_0.2_ppt","species_condition_F_catanatus_0.2_ppt_vs_A_xenica_0.2_ppt")))

#res1 <- results(dds, tidy=TRUE,contrast=list(c("speciesF_catanatus.condition15_ppt")))
#res2 <- results(dds,tidy=TRUE,contrast=list(c("species_condition_L_parva_transfer_vs_A_xenica_0.2_ppt")))

counts_table <- counts(dds, normalized=TRUE )



dim(counts_table)
dim(res_BWvFW)
dim(res_BWvTR)
dim(res_FWvTR)
IDs<-rownames(counts_table)
counts_table<-cbind(counts_table,IDs)
# res_BWvFW
IDs<-res_BWvFW$row
res_BWvFW <- cbind(res_BWvFW,IDs)
counts_stats <- merge(res_BWvFW,counts_table,by="IDs")
dim(counts_stats)
rownames(counts_stats)<-IDs
head(counts_stats)
counts_stats<-counts_stats[-c(2)]
counts_stats<-counts_stats[order(counts_stats$padj),]
head(counts_stats)
names(counts_stats)[names(counts_stats) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(counts_stats)[names(counts_stats) == 'baseMean'] <- 'baseMean-15ppt-v-0.2ppt'
names(counts_stats)[names(counts_stats) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(counts_stats)[names(counts_stats) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(counts_stats)[names(counts_stats) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(counts_stats)[names(counts_stats) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# res_BWvTR
IDs<-res_BWvTR$row
res_BWvTR <- cbind(res_BWvTR,IDs)
counts_stats <- merge(res_BWvTR,counts_stats,by="IDs")
dim(counts_stats)
head(counts_stats)
counts_stats<-counts_stats[-c(2)]
names(counts_stats)[names(counts_stats) == 'padj'] <- 'padj-15ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'baseMean'] <- 'baseMean-15ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'lfcSE'] <- 'lfcSE-15ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'stat'] <- 'stat-15ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'pvalue'] <- 'pvalue-15ppt-v-transfer'
# res_FWvTR
IDs<-res_FWvTR$row
res_FWvTR <- cbind(res_FWvTR,IDs)
counts_stats <- merge(res_FWvTR,counts_stats,by="IDs")
dim(counts_stats)
head(counts_stats)
counts_stats<-counts_stats[-c(2)]
names(counts_stats)[names(counts_stats) == 'padj'] <- 'padj-0.2ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'baseMean'] <- 'baseMean-0.2ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'log2FoldChange'] <- 'log2FoldChange-0.2ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'lfcSE'] <- 'lfcSE-0.2ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'stat'] <- 'stat-0.2ppt-v-transfer'
names(counts_stats)[names(counts_stats) == 'pvalue'] <- 'pvalue-0.2ppt-v-transfer'
counts_stats<-counts_stats[order(counts_stats[,7]),]
head(counts_stats)



write.csv(counts_stats,"~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_results/Ensembl_15killifish_normcounts_3-wayexperimentalconditionstats.csv")

# ============================================
#
# biomart annotation?
# 
# ============================================

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
#ensembl = useDataset("drerio_gene_ensembl",mart=ensembl)
#ensembl = useDataset("amexicanus_gene_ensembl",mart=ensembl)
#ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
ensembl_proteinID = rownames(counts_table)
length(ensembl_proteinID)
query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','go_id','description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
query$entrezgene[query$entrezgene==""] <- "NA"
query_platy<-query
#qery_zlib<-query
#query_cave<-query
#query_mouse<-query
write.csv(query,"~/Documents/UCDavis/Whitehead/reference/biomart_playfish_query.csv")
#ann<-read.csv("~/Documents/UCDavis/Whitehead/counts_annotations.csv")
ann<-counts[,c(2,3,4,5)]
colnames(ann)<-c("gene","scaffold","product","geneID")
#colnames(ann) <- c("x","gene","annotation")



# ============================================
#
# Genes of Interest
# This was very helpful:
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf
# 
# ============================================
# 
#goi <-res$row[c(1:30000)]
#goi <-res$row[c(1:3)]
# Andrew's genes of interest DG/NCBI
# Funhe2EKm029929 XM_012870449.1
# zymogen granule membrane protein 16
goi <- res_BWvTR$row[res_BWvTR$row == "ENSFHEP00000007220.1"]

goi <- res$row[res$row == "XP_012725903.1"]
# Funhe2EKm029931 XM_012870466.1
# zymogen granule membrane protein 16
goi <- res$row[res$row == "XP_012725920.1"]
# solute carrier family 12 member 3-like (removed) 
# Funhe2EKm006896 XM_012852549.1
goi <- res$row[res$row == "XP_012708003.1"]
# chloride channel, voltage-sensitive 2 (clcn2), transcript variant X2 (removed)
# Funhe2EKm024148 XM_012863211.1
goi <- res$row[res$row == "XP_012718665.1"]
# ATP-sensitive inward rectifier potassium channel 1 
# Funhe2EKm001965 XM_012866790.1
goi <- res$row[res$row == "XP_012722244.1"]
# inward rectifier potassium channel 2
#Funhe2EKm023780 XM_012862821.1
goi <- res$row[res$row == "XP_012718275.1"]




# ============================================
# cftr
goi <- res$row[res$row == "XP_012719100.1"]
# polyamine-modulated factor 1-like
goi <- res$row[res$row == "XP_012727384.1"]
# sodium/potassium/calcium exchanger 1 isoform X2
goi <- res$row[res$row == "XP_012706756.1"]
# septin-2B isoform X2
goi <- res$row[res$row == "XP_012716423.1"]
# CLOCK-interacting pacemaker-like
goi <- res$row[res$row == "XP_012722124.1"]
# vasopressin V2 receptor-like
goi <- res$row[res$row == "XP_012721985.1"]
# aquaporin-3 KEEP THIS
goi <- res$row[res$row == "XP_012716807.1"]
# sodium/potassium-transporting ATPase subunit beta-1-interacting protein 1
goi <- res$row[res$row == "XP_012716889.1"]
# septin-2B isoform X2
goi <- res$row[res$row == "XP_012716423.1"]
# otopetrin-1
goi <- res$row[res$row == "XP_012717582.1"]
# claudin-15-like
goi <- res$row[res$row == "XP_012715395.1"]
# claudin-1-like
goi <- res$row[res$row == "XP_012716345.1"]
# claudin 34
goi <- res$row[res$row == "XP_012716562.1"]
# claudin-3-like
goi <- res$row[res$row == "XP_012727928.1"]


#salinity_proteins<-c(207,224,383,420,902,1561,1997,2378,2666,2758,2808,2997,403,492,84,199,469,364,271,217,277,332,532,607,724,848,959,1580,2758,2806,2827,2997,596,672,688,690,778,1130,1266,1955,364,814,1004,1015,2134,614,784,805,1196,1199,1891,1892,1902,2197,2670,2999,3000,248,542,692,1100,1487,2074,2080,2679,861,685,832,1512,1995,249,1333,1478,1957,232,528,280,627,700,2301,373,548,575,1405,383,969,1277,1308,2416,2885,186,1635,189,402,2046,1333,386,186,241,913,1655,596,598,663,672,688,837,1540,2895,257,296,298,387,446,839,1496,231,641)
#length(salinity_proteins)
#goi <-res$row[salinity_proteins]
#goi <- res$row
stopifnot(all(goi %in% names(dds)))
goi
tcounts <- t(log2((counts(dds[goi, ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
  merge(colData(dds), ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(goi)+1):ncol(.))
tcounts_test <- merge(tcounts,ann,by="gene")
tcounts_test %>% select(Row.names, group, species, clade, physiology, condition, gene, product, scaffold,geneID, expression) %>% head %>% knitr::kable()

# --------
# nested facets
# --------

library(gridExtra)

C1<-ggplot(tcounts_test %>%
             filter(clade=='Clade1'),
           aes(condition, expression)) + 
        geom_point(aes(color=physiology)) +
        stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
        facet_grid(~product~species,scales='free_y',labeller=) +
        stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", aes(color=physiology), width=0.2) +
        theme_bw() +
        theme(legend.position="bottom",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x="salinity treatment", 
             y="Expression (log normalized counts)")+
        ggtitle("Clade 1")

#plot(C1)
C2<-ggplot(tcounts_test %>%
             filter(clade=='Clade2'),
           aes(condition, expression)) + 
  geom_point(aes(color=physiology)) +
  stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
  facet_grid(~product~species,scales='free_y',labeller=) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", aes(color=physiology), width=0.2) +
  theme_bw() +
  theme(legend.position="bottom",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="salinity treatment")+
  ggtitle("Clade 2")
#plot(C2)
C3<-ggplot(tcounts_test %>%
             filter(clade=='Clade3'),
           aes(condition, expression)) + 
  geom_point(aes(color=physiology)) +
  stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
  facet_grid(~product~species,scales='free_y',labeller=) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", aes(color=physiology), width=0.2) +
  theme_bw() +
  theme(legend.position="bottom",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="salinity treatment")+
  ggtitle("Clade 3")
#plot(C3)

grid.arrange(C1,C2,C3,ncol=3)

# ----------
# make a pdf
# ----------

pdf("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/multi-ggplot2-catalog_salinity_31Oct2018.pdf")
for (i in goi) {
  p <- ggplot(filter(tcounts_test, gene==i), aes(condition, expression, fill=physiology)) +  geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~product~clade+species,scales='free_y') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    labs(x="salinity treatment")
  print(p)
}
dev.off()



plotPCAWithSampleNames(log_cds,intgroup="clade",ntop=40000)

counts_table = counts(dds, normalized=TRUE )

dim(counts_table)






filtered_norm_counts<-subset(counts_table,rownames(counts_table) %in% goi)
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
dim(filtered_norm_counts)
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)

res.1<-results(dds,contrast=c("condition","15_ppt","0.2_ppt"))
dim(res.1)
res.2<-results(dds,contrast=c("clade","Clade1","Clade3"))

res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
res1_filtered <- subset(res1_ordered,rownames(res1_ordered) %in% goi)
#res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
#res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
goi
res1_filtered<-cbind(res1_filtered,id)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
#res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="Salinity exposure (15 ppt vs. 0.2 ppt) (padj<0.05))",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(rownames(res.1) %in% goi,"blue","gray67"),
     main="Salinity exposure (15 ppt vs. 0.2 ppt) (known salinity responsive genes = blue))",xlim=c(1,15),pch=20,cex=1)
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
heatmap.2(d, main="16 species killifish, known salinity-responsive genes",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)