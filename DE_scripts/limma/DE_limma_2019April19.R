library(DESeq2)
library(limma)
library('edgeR')
library("RColorBrewer")
library(gplots)
library(lattice)
library("vsn")
#setwd("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_scripts/limma")

# This file contains filtered counts with Experimental Design Info in the last 5 rows
if(!file.exists('../../../../mean1_species_counts_design.csv')){
  download.file("https://osf.io/j9ewk/download",'mean1_species_counts_design.csv')
}


#counts_design <- read.csv("../../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = TRUE)
#counts_design <- read.csv("../../../nonzero_clade_physiology_counts_design.csv",stringsAsFactors = TRUE)
#counts_design <- read.csv("../../../nozero_species_counts_design.csv",stringsAsFactors = TRUE)
counts_design <- read.csv("../../../mean1_species_counts_design.csv",stringsAsFactors = TRUE)
dim(counts_design)
#[1] 9111  129

# -----------------------
# Format design and counts matrix
# Drop columns with no data
# -----------------------

design <- counts_design[counts_design$Ensembl == 'Empty',]
#design$type <- c("species","native_salinity","clade","group","condition")
drops <- c("X","Ensembl",
           "F_zebrinus_BW_1.quant","F_zebrinus_BW_2.quant",
           "F_zebrinus_FW_1.quant","F_zebrinus_FW_2.quant",
           "F_notti_FW_1.quant","F_notti_FW_2.quant",
           "F_sciadicus_BW_1.quant","F_sciadicus_FW_1.quant","F_sciadicus_FW_2.quant")
transfer_drops <- c("F_sciadicus_transfer_1.quant","F_rathbuni_transfer_1.quant","F_rathbuni_transfer_2.quant","F_rathbuni_transfer_3.quant",
                    "F_grandis_transfer_1.quant","F_grandis_transfer_2.quant","F_grandis_transfer_3.quant",
                    "F_notatus_transfer_1.quant","F_notatus_transfer_2.quant","F_notatus_transfer_3.quant",
                    "F_parvapinis_transfer_1.quant","F_parvapinis_transfer_2.quant",
                    "L_goodei_transfer_1.quant","L_goodei_transfer_2.quant","L_goodei_transfer_3.quant",
                    "F_olivaceous_transfer_1.quant","F_olivaceous_transfer_2.quant",
                    "L_parva_transfer_1.quant","L_parva_transfer_2.quant","L_parva_transfer_3.quant",
                    "F_heteroclitusMDPP_transfer_1.quant","F_heteroclitusMDPP_transfer_2.quant","F_heteroclitusMDPP_transfer_3.quant",
                    "F_similis_transfer_1.quant","F_similis_transfer_2.quant","F_similis_transfer_3.quant",
                    "F_diaphanus_transfer_1.quant","F_diaphanus_transfer_2.quant",
                    "F_chrysotus_transfer_1.quant","F_chrysotus_transfer_2.quant",
                    "A_xenica_transfer_1.quant","A_xenica_transfer_2.quant","A_xenica_transfer_3.quant" ,
                    "F_catanatus_transfer_1.quant","F_catanatus_transfer_2.quant",
                    "F_heteroclitusMDPL_transfer_1.quant","F_heteroclitusMDPL_transfer_2.quant","F_heteroclitusMDPL_transfer_3.quant")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
design <- design[ , !(names(design) %in% transfer_drops)]
counts <- counts[ , !(names(counts) %in% transfer_drops)]
dim(design)
#[1]  5 81
dim(counts)
#[1] 9106   81

# --------------------
# design cateogories
# --------------------

species<-as.character(unlist(design[1,]))
physiology<-as.character(unlist(design[2,]))
clade<-as.character(unlist(design[3,]))
condition<-as.character(unlist(design[5,]))
condition_physiology<-as.vector(paste(condition,physiology,sep="."))
cols<-colnames(counts)
ExpDesign <- data.frame(row.names=cols,
                        condition=condition,
                        physiology = physiology,
                        clade = clade,
                        species = species,
                        sample=cols)
ExpDesign
design = model.matrix( ~ physiology + condition + physiology:condition, ExpDesign)

colnames(design)
# check rank of matrix
Matrix::rankMatrix( design )
dim(design)

# ---------------

# stabilize variance

# ---------------

counts_round <- round(data.matrix(counts), digits=0)
#counts_round <- head(counts_round, n = 19000)
dim(counts_round)
plot(colSums(t(counts_round)) )
boxplot(counts_round, las = 2, main = "raw counts")

# transform with DESeq to stabilize variance

dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = design)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds,fitType='local')
vsd <- varianceStabilizingTransformation(dds,fitType='local')

plot(rank(rowMeans(counts(dds))),
     genefilter::rowVars(log2(counts(dds)+1)), main="log2(x+1) transform")
plot(rank(rowMeans(assay(vsd))), genefilter::rowVars(assay(vsd)),
     main="VST")

ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

plot(colSums(t(assay(vsd))) )

# this takes a long time
#rld <- rlog(dds, blind=FALSE)
#meanSdPlot(assay(rld))
# ---------------

# DE analysis

# ---------------
tf <- assay(vsd) 
tf[tf < 0] <- NA
tf[is.na(tf)] <- 0

#genes = DGEList(count = assay(vsd), group = condition_physiology)
# Error: Negative counts not allowed
#genes = DGEList(count = tf, group = condition_physiology)
genes = DGEList(count = counts_round, group = condition_physiology)
#genes = calcNormFactors( genes )
vobj = voom( genes, design, plot=TRUE)
lcpm <- cpm(genes$counts, log = TRUE)
boxplot(lcpm, las = 2, main = "After limma-voom Normalization")



col.group <- condition_physiology
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)


pca = prcomp(t(lcpm))
names = colnames(lcpm)
fac= factor(condition_physiology)
colours = c("red","blue","green","orange")
xyplot(
  PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
  panel=function(x, y, ...) {
    panel.xyplot(x, y, ...);
    ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=1)
  },
  aspect = "fill", col=colours
  #main = draw.key(key = list(rect = list(col = list(col=colours), text = list(levels(fac)), rep = FALSE)))
)


#vobj = voom( genes, design, plot=TRUE)
vwts <- voomWithQualityWeights(genes, design=design, normalization="quantile", plot=TRUE)

corfit <- duplicateCorrelation(vobj,design,block=ExpDesign$clade)

corfit$consensus
# [1] 0.01488349

fitRan <- lmFit(vobj,design,block=ExpDesign$clade,correlation=corfit$consensus)
fitRan <- eBayes(fitRan)
topTable(fitRan,number=50,coef=ncol(design))
