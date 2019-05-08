library(limma)
library('edgeR')
library("RColorBrewer")
library(gplots)
library(lattice)
library("vsn")
setwd("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_scripts/limma")

# This file contains filtered counts with Experimental Design Info in the last 5 rows
if(!file.exists('../../../../Ensembl_species_counts_designfactors.csv')){
  download.file("https://osf.io/7vp38/download",'Ensembl_species_counts_designfactors.csv')
}

counts_design <- read.csv("../../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)

dim(counts_design)
#[1] 30471   130

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
gene.names<-rownames(counts)

dim(counts)
class(counts)
design[] <- lapply( design, factor)

# --------------------
# design cateogories
# --------------------

species<-as.character(unlist(design[1,]))
physiology<-as.character(unlist(design[2,]))
clade<-as.character(unlist(design[3,]))
condition<-as.character(unlist(design[5,]))
condition_physiology<-as.vector(paste(condition,physiology,sep="."))
condition_physiology_clade <- as.vector(paste(condition_physiology,clade,sep="."))
condition_physiology_clade <- as.vector(paste("group",condition_physiology_clade,sep=""))
cols<-colnames(counts)
ExpDesign <- data.frame(row.names=cols,
                        condition=condition,
                        physiology = physiology,
                        clade = clade,
                        species = species,
                        sample=cols)
ExpDesign
design = model.matrix( ~0 + physiology*condition*clade, ExpDesign)
colnames(design)
# check rank of matrix
Matrix::rankMatrix( design )
dim(design)
design2<-model.matrix(~0 + genes$samples$group)
colnames(design2)<- gsub("genes\\$samples\\$","",colnames(design2))
colnames(design2)
rownames(design2)<-rownames(design)
# ---------------

# DE analysis

# ---------------

counts<-as.matrix(as.data.frame(sapply(counts, as.numeric)))
rownames(counts)<-gene.names
class(counts)

keep<-filterByExpr(counts,design = design,group=condition_physiology,min.count = 10, min.total.count = 100)
counts.filt <- counts[keep,]
dim(counts.filt)
write.csv(counts.filt,"../../../21k_counts_filt_30April2019.csv")

genes = DGEList(count = counts.filt, group = condition_physiology_clade)
#genes = calcNormFactors( genes )


vobj = voom( genes, design, plot=TRUE)
lcpm <- cpm(genes$counts, log = TRUE)
boxplot(lcpm, las = 2, main = "")
plot(colSums(t(lcpm)))

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

vwts <- voomWithQualityWeights(genes, design=design, normalization="quantile", plot=TRUE)

corfit <- duplicateCorrelation(vobj,design,block=ExpDesign$species)

corfit$consensus
#[1] 0.751381

fitRan <- lmFit(vobj,design,block=ExpDesign$species,correlation=corfit$consensus)
colnames(coef(fitRan))
# [1] "(Intercept)"                             "physiologyM"                            
# [3] "condition15_ppt"                         "cladeClade2"                            
# [5] "cladeClade3"                             "physiologyM:condition15_ppt"            
# [7] "physiologyM:cladeClade2"                 "physiologyM:cladeClade3"                
# [9] "condition15_ppt:cladeClade2"             "condition15_ppt:cladeClade3"            
# [11] "physiologyM:condition15_ppt:cladeClade2" "physiologyM:condition15_ppt:cladeClade3"

cont.matrix <- makeContrasts(condition15_ppt - ,levels=make.names(ExpDesign))

cont.matrix
tmp <- contrasts.fit(fitRan, cont.matrix)
fitRan <- eBayes(tmp)

stats<-topTable(fitRan,number=50,sort.by="p")
topTable(fitRan,number=50,sort.by="p")

