library(limma)
library(ggplot2)
library('edgeR')
library("RColorBrewer")
library(gplots)
library(biomaRt)
library(lattice)
library("vsn")
library(dplyr)
library(knitr)
library(kableExtra)
library(pheatmap)
library("SummarizedExperiment")
library("emmeans")
setwd("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_scripts/limma")
dir<-("~/Documents/UCDavis/Whitehead/")
# This file contains filtered counts with Experimental Design Info in the last 5 rows
if(!file.exists('~/Documents/UCDavis/Whitehead/kfish_expression_July2019/Ensembl_species_counts_designfactors.csv')){
  download.file("https://osf.io/7vp38/download",'Ensembl_species_counts_designfactors.csv')
}

counts_design <- read.csv('~/Documents/UCDavis/Whitehead/kfish_expression_July2019/Ensembl_species_counts_designfactors.csv',stringsAsFactors = FALSE)

dim(counts_design)
#[1] 31595   130

# -----------------------
# Format design and counts matrix
# Drop columns with no data
# -----------------------

design <- counts_design[counts_design$Ensembl == 'Empty',]
#design$type <- c("species","native_salinity","clade","group","condition")
drops <- c("X","Ensembl",
           "F_zebrinus_BW_1.quant","F_zebrinus_BW_2.quant",
           "F_zebrinus_FW_1.quant","F_zebrinus_FW_2.quant",
           "F_nottii_FW_1.quant","F_nottii_FW_2.quant",
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
# design categories
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
# used for pairwise contrasts
#form<-as.formula("~0 + physiology*condition*clade")
form<-as.formula("~physiology*condition*clade")
design = model.matrix(form, ExpDesign)
#group <- interaction(physiology, condition, clade)
#mm <- model.matrix(~0 + group)
colnames(design)
# check rank of matrix
Matrix::rankMatrix( design )
dim(design)
clade <- ExpDesign$clade
physiology <- ExpDesign$physiology

# ============================================
# biomart annotation
# https://uswest.ensembl.org/Fundulus_heteroclitus/Info/Index
# ============================================

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
# test
# id <- c("ENSGMOP00000000001.1","ENSGMOP00000000002.1","ENSGMOP00000000003.1")
ensembl_proteinID = rownames(counts)
length(ensembl_proteinID)
# can take >5 min
# do only once
ann<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','description'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
head(ann)
dim(ann)
length(unique(ann$ensembl_peptide_id))
ann <- ann[!duplicated(ann[,c(1)]),]
dim(ann)

# ---------------

# DE analysis

# ---------------

counts<-as.matrix(as.data.frame(sapply(counts, as.numeric)))
rownames(counts)<-gene.names
class(counts)
#test<-counts %>% drop_na()
#test<-as.matrix(test)
lcom_unfilt<-log2(counts+1)
plot(colSums(t(lcom_unfilt)))

keep<-filterByExpr(counts,design = design,group=condition_physiology,min.count = 10, min.total.count = 100)
counts.filt <- counts[keep,]
dim(counts.filt)
#write.table(counts.filt,"~/Documents/UCDavis/Whitehead/exp.tsv",sep="\t",quote=FALSE)

norm_counts <- read.csv("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/normalized_counts.csv",stringsAsFactors = FALSE, header = TRUE,row.names = "Ensembl")
dim(norm_counts)

# ============================================
#
# Genes of Interest
# This was very helpful:
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf
# 
# ============================================

# ---------------------------
# Andrew Whitehead's genes of interest:
# ---------------------------

# zymogen granule membrane protein 16
# Funhe2EKm029929
countsfilt <- counts.filt
countsfilt$row <- rownames(countsfilt)
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000007220.1"]
goi
# zymogen granule membrane protein 16
# Funhe2EKm029931
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000025841"]
goi
# solute carrier family 12 member 3-like (removed) 
# Funhe2EKm006896
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000009214"]
goi
# chloride channel, voltage-sensitive 2 (clcn2), transcript variant X2 (removed)
# Funhe2EKm024148
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000019510"]
goi
# ATP-sensitive inward rectifier potassium channel 1 
# Funhe2EKm001965
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000015383"]
goi
# inward rectifier potassium channel 2
#Funhe2EKm023780
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000009753"]

# --------------------------------
# other salinity genes of interest
# --------------------------------

# ============================================
# aquaporin-3
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000006725"]
goi
# cftr
# Funhe2EKm024501
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000008393"]
goi
# polyamine-modulated factor 1-like
# Funhe2EKm031049
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000013324"]
goi
# sodium/potassium/calcium exchanger 2
# ENSFHEP00000034177
# Funhe2EKm025497
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000034177"]
goi
# septin-2B isoform X2
# ENSFHEP00000015765
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000015765"]
goi
# CLOCK-interacting pacemaker-like
# ENSFHEP00000017303
# Funhe2EKm026846
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000017303"]
goi
# vasopressin V2 receptor-like
# Funhe2EKm026721
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000000036"]
goi
# sodium/potassium-transporting ATPase subunit beta-1-interacting protein 1
# ENSFHEP00000031108
# Funhe2EKm001174
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000031108"]
goi
# ENSFHEP00000023396
# sodium/potassium ATPase alpha subunit isoform 2 
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000023396"]
goi
# septin-2
# Funhe2EKm012182
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000016853"]
goi
# otopetrin-2
# Funhe2EKm035427
goi <- countsfiltW$row[countsfilt$row == "ENSFHEP00000026411"]
goi
# claudin 8
# Funhe2EKm037718
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000006282"]
goi
# claudin 4
# ENSFHEP00000003908
# Funhe2EKm007149
goi <- countsfilt$row[countsfilt$row == "ENSFHEP00000003908"]
goi
all_goi<-c("ENSFHEP00000023396","ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
           "ENSFHEP00000015383","ENSFHEP00000009753","ENSFHEP00000006725","ENSFHEP00000008393",
           "ENSFHEP00000013324","ENSFHEP00000001609","ENSFHEP00000013324","ENSFHEP00000034177",
           "ENSFHEP00000015765","ENSFHEP00000017303","ENSFHEP00000000036","ENSFHEP00000031108",
           "ENSFHEP00000016853","ENSFHEP00000003908")

# -------------------
# proceed with DE
# -------------------
# log counts before DE
boxplot(log(counts.filt+1), las = 2, main = "")

#write.table(counts.filt,"~/Documents/UCDavis/Whitehead/exp.tsv",sep = "\t", quote = F, row.names = F)

genes = DGEList(count = counts.filt, group = condition_physiology_clade)
genes = calcNormFactors(genes)

# write normalized counts
dir <- "~/Documents/UCDavis/Whitehead/"
tmp <- as.data.frame(cpm(genes))
tmp$Ensembl <- rownames(tmp)
tmp <- dplyr::select(tmp, Ensembl, everything())
#write.csv(tmp, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/normalized_counts.csv"), quote = F, row.names = F)

vobj = voom( genes, design, plot=TRUE)
lcpm <- cpm(genes$counts, log = TRUE)

# log counts after DE

boxplot(lcpm, las = 2, main = "")
plot(colSums(t(lcpm)))

vwts <- voomWithQualityWeights(genes, design=design, normalization="quantile", plot=TRUE)
corfit <- duplicateCorrelation(vobj,design,block=ExpDesign$species)
corfit$consensus
#[1] 0.758966
fitRan <- lmFit(vobj,design,block=ExpDesign$species,correlation=corfit$consensus)

# -----------------
# PCA, all counts
# -----------------

# This is the PCA with all counts, not filtered. The dimensions of the counts table are listed below. Row=genes, Columns=samples
# normal full counts
x <- data.matrix(cpm(genes))
dim(x)
x <- x+1
log_x<-log(x)
names<-colnames(log_x)
pca = prcomp(t(log_x))
summary(pca)
fac = factor(physiology)
colours = function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])}
#mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(pca$x[,1:2], 
     col=colours(clade), 
     pch = c(16, 2, 9)[as.numeric(as.factor(physiology))],
     cex=2,
     xlab="PC1",
     ylab="PC2",
     cex.lab=2,
     cex.axis = 2)
#legend(140,100,legend=c("Clade 1","Clade 2","Clade 3"),col=rainbow(length(unique(clade))),cex=0.75, pch=19)
#legend(140,-67,legend=c("Freshwater","Marine"),cex=0.75,pch=c(16, 2, 9))
legend(-75,50,legend=c("Clade 1","Clade 2","Clade 3"),col=rainbow(length(unique(clade))),cex=0.75, pch=19)
legend(-75,25,legend=c("Freshwater","Marine"),cex=0.75,pch=c(16, 2, 9))

# -------------------
# coefficient names
# -------------------
colnames(coef(fitRan))
dir <- "~/Documents/UCDavis/Whitehead/"
#> colnames(coef(fitRan))
#[1] "(Intercept)"                             "physiologyM"                            
#[3] "condition15_ppt"                         "cladeClade2"                            
#[5] "cladeClade3"                             "physiologyM:condition15_ppt"            
#[7] "physiologyM:cladeClade2"                 "physiologyM:cladeClade3"                
#[9] "condition15_ppt:cladeClade2"             "condition15_ppt:cladeClade3"            
#[11] "physiologyM:condition15_ppt:cladeClade2" "physiologyM:condition15_ppt:cladeClade3"

# ---------------------
# salinity (condition) main effect
# ---------------------
contrasts_salinity_main <- makeContrasts("condition15_ppt", "physiologyM.condition15_ppt","condition15_ppt.cladeClade2","condition15_ppt.cladeClade3","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_salinity_main) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_salinity_main)
vfit <- eBayes(vfit)
main_salinity<-topTableF(vfit,n = Inf, sort.by = "F")
sig_main_salinity <- main_salinity[main_salinity$adj.P.Val<0.01,]
dim(sig_main_salinity)
dim(main_salinity)
ann_salinity <- merge(main_salinity,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_salinity <- ann_salinity[!is.na(ann_salinity$adj.P.Val),]
dim(ann_salinity)
ann_salinity <- ann_salinity[order(ann_salinity$adj.P.Val,decreasing = FALSE), ]
rownames(ann_salinity) <- ann_salinity$Row.names
ann_salinity <- ann_salinity[,-1]
#write.csv(ann_salinity, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/main_salinity.csv"), quote = F, row.names = T)

# ---------------------
# physiology main effect
# ---------------------
contrasts_physiology_main <- makeContrasts("physiologyM", "physiologyM.cladeClade2","physiologyM.cladeClade3","physiologyM.condition15_ppt","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_physiology_main) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_physiology_main)
vfit <- eBayes(vfit)
main_physiology<-topTableF(vfit,n = Inf, sort.by = "F")
sig_main_physiology <- main_physiology[main_physiology$adj.P.Val<0.05,]
dim(sig_main_physiology)
dim(main_physiology)
ann_physiology <- merge(main_physiology,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_physiology <- ann_physiology[!is.na(ann_physiology$adj.P.Val),]
dim(ann_physiology)
ann_physiology <- ann_physiology[order(ann_physiology$adj.P.Val,decreasing = FALSE), ]
rownames(ann_physiology) <- ann_physiology$Row.names
ann_physiology <- ann_physiology[,-1]
#write.csv(ann_physiology, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/main_physiology.csv"), quote = F, row.names = T)


# ---------------------
# clade main effect
# ---------------------
contrasts_clade_main <- makeContrasts("cladeClade2","cladeClade3","physiologyM.cladeClade2","physiologyM.cladeClade3","condition15_ppt.cladeClade2","condition15_ppt.cladeClade3","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_clade_main) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_clade_main)
vfit <- eBayes(vfit)
main_clade<-topTableF(vfit,n = Inf, sort.by = "F")
sig_main_clade <- main_clade[main_clade$adj.P.Val<0.05,]
dim(sig_main_clade)
dim(main_clade)
ann_clade <- merge(main_clade,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_clade <- ann_clade[!is.na(ann_clade$adj.P.Val),]
dim(ann_clade)
ann_clade <- ann_clade[order(ann_clade$adj.P.Val,decreasing = FALSE), ]
rownames(ann_clade) <- ann_clade$Row.names
ann_clade <- ann_clade[,-1]
#write.csv(ann_clade, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/main_clade.csv"), quote = F, row.names = T)


# ---------------------
# three-way interaction
# ---------------------
contrasts_threeway <- makeContrasts("physiologyM.condition15_ppt.cladeClade2", "physiologyM.condition15_ppt.cladeClade3", levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_threeway) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_threeway)
vfit <- eBayes(vfit)
main_threeway<-topTableF(vfit,n = Inf, sort.by = "F")
sig_threeway <- main_threeway[main_threeway$adj.P.Val<0.05,]
dim(sig_threeway)
dim(main_threeway)
ann_threeway <- merge(main_threeway,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_threeway <- ann_threeway[!is.na(ann_threeway$adj.P.Val),]
dim(ann_threeway)
ann_threeway <- ann_threeway[order(ann_threeway$adj.P.Val,decreasing = FALSE), ]
rownames(ann_threeway) <- ann_threeway$Row.names
ann_threeway <- ann_threeway[,-1]
#write.csv(ann_threeway, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/threeway.csv"), quote = F, row.names = T)

# ---------------------
# salinity x clade two-way interaction
# ---------------------
contrasts_salinity_clade_interaction <- makeContrasts("condition15_ppt.cladeClade2","condition15_ppt.cladeClade3","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_salinity_clade_interaction) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_salinity_clade_interaction)
vfit <- eBayes(vfit)
main_salinity_clade_interaction<-topTableF(vfit,n = Inf, sort.by = "F")
sig_salinity_clade_interaction <- main_salinity_clade_interaction[main_salinity_clade_interaction$adj.P.Val<0.05,]
dim(sig_salinity_clade_interaction)
dim(main_salinity_clade_interaction)
ann_salinity_clade_interaction <- merge(main_salinity_clade_interaction,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_salinity_clade_interaction <- ann_salinity_clade_interaction[!is.na(ann_salinity_clade_interaction$adj.P.Val),]
dim(ann_salinity_clade_interaction)
ann_salinity_clade_interaction <- ann_salinity_clade_interaction[order(ann_salinity_clade_interaction$adj.P.Val,decreasing = FALSE), ]
rownames(ann_salinity_clade_interaction) <- ann_salinity_clade_interaction$Row.names
ann_salinity_clade_interaction <- ann_salinity_clade_interaction[,-1]
#write.csv(ann_salinity_clade_interaction, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/salinity_clade_interaction.csv"), quote = F, row.names = T)

# ---------------------
# salinity x physiology two-way interaction
# ---------------------
contrasts_salinity_physiology_interaction <- makeContrasts("physiologyM.condition15_ppt","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))
rownames(contrasts_salinity_physiology_interaction) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_salinity_physiology_interaction)
vfit <- eBayes(vfit)
main_salinity_physiology_interaction<-topTableF(vfit,n = Inf, sort.by = "F")
sig_salinity_physiology_interaction <- main_salinity_physiology_interaction[main_salinity_physiology_interaction$adj.P.Val<0.05,]
dim(sig_salinity_physiology_interaction)
dim(main_salinity_physiology_interaction)
ann_salinity_physiology_interaction <- merge(main_salinity_physiology_interaction,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_salinity_physiology_interaction <- ann_salinity_physiology_interaction[!is.na(ann_salinity_physiology_interaction$adj.P.Val),]
dim(ann_salinity_physiology_interaction)
ann_salinity_physiology_interaction <- ann_salinity_physiology_interaction[order(ann_salinity_physiology_interaction$adj.P.Val,decreasing = FALSE), ]
rownames(ann_salinity_physiology_interaction) <- ann_salinity_physiology_interaction$Row.names
ann_salinity_physiology_interaction <- ann_salinity_physiology_interaction[,-1]
#write.csv(ann_salinity_physiology_interaction, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/salinity_physiology_interaction.csv"), quote = F, row.names = T)

# ---------------------
# physiology x clade two-way interaction
# ---------------------
contrasts_clade_physiology_interaction <- makeContrasts("physiologyM.cladeClade2","physiologyM.cladeClade3","physiologyM.condition15_ppt.cladeClade2","physiologyM.condition15_ppt.cladeClade3",levels=make.names(colnames(coef(fitRan))))

rownames(contrasts_clade_physiology_interaction) <- colnames(coef(fitRan))
vfit <- contrasts.fit(fitRan, contrasts = contrasts_clade_physiology_interaction)
vfit <- eBayes(vfit)
main_clade_physiology_interaction<-topTableF(vfit,n = Inf, sort.by = "F")
sig_clade_physiology_interaction <- main_clade_physiology_interaction[main_clade_physiology_interaction$adj.P.Val<0.05,]
dim(sig_clade_physiology_interaction)
dim(main_clade_physiology_interaction)
ann_clade_physiology_interaction <- merge(main_clade_physiology_interaction,ann,all=TRUE,by.x = "row.names", by.y  = "ensembl_peptide_id")
ann_clade_physiology_interaction <- ann_clade_physiology_interaction[!is.na(ann_clade_physiology_interaction$adj.P.Val),]
dim(ann_clade_physiology_interaction)
ann_clade_physiology_interaction <- ann_clade_physiology_interaction[order(ann_clade_physiology_interaction$adj.P.Val,decreasing = FALSE), ]
rownames(ann_clade_physiology_interaction) <- ann_clade_physiology_interaction$Row.names
ann_clade_physiology_interaction <- ann_clade_physiology_interaction[,-1]
#write.csv(ann_clade_physiology_interaction, file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/clade_physiology_interaction.csv"), quote = F, row.names = T)



#======================

#-------------------------------------
# significant genes
#-------------------------------------

dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)
dim(sig_threeway)
dim(sig_main_clade)
dim(sig_main_physiology)
dim(sig_main_salinity)

# make a summary table for Andrew:
# ---------------------------------
# cols:
# Gene ID
# gene name
# gene descriptions
# adj.P.Val.Clade
# adj.P.Val.Salinity
# adj.P.Val.Physiology
# adj.P.Val.Salinity.Physiology
# adj.P.Val.Physiology.Clade
# adj.P.Val.Salinity.Clade
# adj.P.Val.3way
# norm counts
# ---------------------------------
new <- merge(ann_clade,ann_physiology,by= "row.names")
colnames(new)[13] <- "adj.P.Val.Clade"
colnames(new)[28] <- "adj.P.Val.Physiology"
new <- new[ -c(2:12, 19:27,29:33) ]
colnames(new)[1] <- "ID"
new <- merge(new,ann_salinity,by.x= "ID",by.y="row.names")
colnames(new)[18] <- "adj.P.Val.Salinity"
new <- new[ -c(9:17,19:23) ]
new <- new[,c(1,3:7,2,8:9)]
new <- merge(new,ann_clade_physiology_interaction,by.x= "ID",by.y="row.names")
colnames(new)[17] <- "adj.P.Val.Physiology.Clade"
new <- new[ -c(10:16,18:22) ]
new <- merge(new,ann_salinity_physiology_interaction,by.x= "ID",by.y="row.names")
colnames(new)[17] <- "adj.P.Val.Salinity.Physiology"
new <- new[ -c(11:16,18:22) ]
new <- merge(new,ann_salinity_clade_interaction,by.x= "ID",by.y="row.names")
colnames(new)[19] <- "adj.P.Val.Salinity.Clade"
new <- new[ -c(12:18,20:24) ]
new <- merge(new,ann_threeway,by.x= "ID",by.y="row.names")
colnames(new)[18] <- "adj.P.Val.3way"
new <- new[ -c(13:17,19:23) ]
colnames(new)[2] <- "ensembl_transcript_id"
colnames(new)[3] <- "ensembl_gene_id"
colnames(new)[4] <- "gene_biotype"
colnames(new)[5] <- "external_gene_name"
colnames(new)[6] <- "description"
new <- merge(new,norm_counts,by.x= "ID",by.y="row.names")
# ---------------------------------
# mean expression each species BW
# mean expression each species FW
# x 14
new$F.rath.BW.mean <- log2(rowMeans(new[c(14:16)], na.rm=TRUE)+1)
new$F.rath.FW.mean <- log2(rowMeans(new[c(17:19)], na.rm=TRUE)+1)
new$F.grandis.BW.mean <- log2(rowMeans(new[c(20:22)], na.rm=TRUE)+1)
new$F.grandis.FW.mean <- log2(rowMeans(new[c(23:25)], na.rm=TRUE)+1)
new$F.notatus.BW.mean <- log2(rowMeans(new[c(26:28)], na.rm=TRUE)+1)
new$F.notatus.FW.mean <- log2(rowMeans(new[c(29:31)], na.rm=TRUE)+1)
new$F.parv.BW.mean <- log2(rowMeans(new[c(32:34)], na.rm=TRUE)+1)
new$F.parv.FW.mean <- log2(rowMeans(new[c(35:37)], na.rm=TRUE)+1)
new$L.good.BW.mean <- log2(rowMeans(new[c(38:40)], na.rm=TRUE)+1)
new$L.good.FW.mean <- log2(rowMeans(new[c(41:43)], na.rm=TRUE)+1)
new$F.oli.BW.mean <- log2(rowMeans(new[c(44:46)], na.rm=TRUE)+1)
new$F.oli.FW.mean <- log2(rowMeans(new[c(47:49)], na.rm=TRUE)+1)
new$L.parv.BW.mean <- log2(rowMeans(new[c(50:52)], na.rm=TRUE)+1)
new$L.parv.FW.mean <- log2(rowMeans(new[c(53:55)], na.rm=TRUE)+1)
new$F.hetPP.BW.mean <- log2(rowMeans(new[c(56:58)], na.rm=TRUE)+1)
new$F.hetPP.FW.mean <- log2(rowMeans(new[c(59:61)], na.rm=TRUE)+1)
new$F.sim.BW.mean <- log2(rowMeans(new[c(62:64)], na.rm=TRUE)+1)
new$F.sim.FW.mean <- log2(rowMeans(new[c(65:67)], na.rm=TRUE)+1)
new$F.dia.BW.mean <- log2(rowMeans(new[c(68:69)], na.rm=TRUE)+1)
new$F.dia.FW.mean <- log2(rowMeans(new[c(70:71)], na.rm=TRUE)+1)
new$F.chry.BW.mean <- log2(rowMeans(new[c(72:74)], na.rm=TRUE)+1)
new$F.chry.FW.mean <- log2(rowMeans(new[c(72:77)], na.rm=TRUE)+1)
new$A.xen.BW.mean <- log2(rowMeans(new[c(78:80)], na.rm=TRUE)+1)
new$A.xen.FW.mean <- log2(rowMeans(new[c(81:83)], na.rm=TRUE)+1)
new$F.cat.BW.mean <- log2(rowMeans(new[c(84:86)], na.rm=TRUE)+1)
new$F.cat.FW.mean <- log2(rowMeans(new[c(87:88)], na.rm=TRUE)+1)
new$F.hetPL.BW.mean <- log2(rowMeans(new[c(89:91)], na.rm=TRUE)+1)
new$F.hetPL.FW.mean <- log2(rowMeans(new[c(92:94)], na.rm=TRUE)+1)
# ---------------------------------
# log2FC FW=0 - each species - zero normalized
# log2FC BW=15 - each species

# subtract FW from itself = 0
# subtract FW from BW = 15

new$F.rath.0 <- (new$F.rath.FW.mean - new$F.rath.FW.mean) 
new$F.rath.15 <- (new$F.rath.BW.mean - new$F.rath.FW.mean)
new$F.grandis.0 <- (new$F.grandis.FW.mean - new$F.grandis.FW.mean)
new$F.grandis.15 <- (new$F.grandis.BW.mean - new$F.grandis.FW.mean)
new$F.notatus.0 <- (new$F.notatus.FW.mean - new$F.notatus.FW.mean)
new$F.notatus.15 <- (new$F.notatus.BW.mean - new$F.notatus.FW.mean)
new$F.parv.0 <- (new$F.parv.FW.mean - new$F.parv.FW.mean)
new$F.parv.15 <- (new$F.parv.BW.mean - new$F.parv.FW.mean)
new$L.good.0 <- (new$L.good.FW.mean - new$L.good.FW.mean)
new$L.good.15 <- (new$L.good.BW.mean - new$L.good.FW.mean)
new$F.oli.0 <- (new$F.oli.FW.mean - new$F.oli.FW.mean)  
new$F.oli.15 <- (new$F.oli.BW.mean - new$F.oli.FW.mean)  
new$L.parv.0 <- (new$L.parv.FW.mean - new$L.parv.FW.mean) 
new$L.parv.15 <- (new$L.parv.BW.mean - new$L.parv.FW.mean)
new$F.hetPP.0 <- (new$F.hetPP.FW.mean - new$F.hetPP.FW.mean) 
new$F.hetPP.15 <- (new$F.hetPP.BW.mean - new$F.hetPP.FW.mean)
new$F.sim.0 <- (new$F.sim.FW.mean - new$F.sim.FW.mean)
new$F.sim.15 <- (new$F.sim.BW.mean - new$F.sim.FW.mean)
new$F.dia.0 <- (new$F.dia.FW.mean - new$F.dia.FW.mean)
new$F.dia.15 <- (new$F.dia.BW.mean - new$F.dia.FW.mean)
new$F.chry.0 <- (new$F.chry.FW.mean - new$F.chry.FW.mean)
new$F.chry.15 <- (new$F.chry.BW.mean - new$F.chry.FW.mean)
new$A.xen.0 <- (new$A.xen.FW.mean - new$A.xen.FW.mean)
new$A.xen.15 <- (new$A.xen.BW.mean - new$A.xen.FW.mean)
new$F.cat.0 <- (new$F.cat.FW.mean - new$F.cat.FW.mean)
new$F.cat.15 <- (new$F.cat.BW.mean - new$F.cat.FW.mean)
new$F.hetPL.0 <- (new$F.hetPL.FW.mean - new$F.hetPL.FW.mean)
new$F.hetPL.15 <- (new$F.hetPL.BW.mean - new$F.hetPL.FW.mean)
#write.csv(new,file = file.path("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/kfishosmotictolerance_Summary.Table.csv"), quote = F, row.names = T)

#-------------------
# Set sample order
#-------------------

sample <- c("F.grandis.0","F.grandis.15","F.hetPL.0","F.hetPL.15","F.hetPP.0","F.hetPP.15","F.dia.0","F.dia.15","F.cat.0","F.cat.15","F.rath.0","F.rath.15","F.parv.0","F.parv.15","L.parv.0","L.parv.15","L.good.0","L.good.15","A.xen.0","A.xen.15","F.chry.0","F.chry.15","F.sim.0","F.sim.15","F.notatus.0","F.notatus.15","F.oli.0","F.oli.15")
clade <- c("clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade2","clade2","clade2","clade2","clade2","clade2","clade3","clade3","clade3","clade3","clade3","clade3","clade3","clade3","clade3","clade3")
physiology <- c("M","M","M","M","M","M","M","M","F","F","F","F","M","M","M","M","F","F","M","M","M","M","M","M","F","F","F","F")
salinity <- c("0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt")
sample_label_df <- data.frame(salinity,physiology,clade)
rownames(sample_label_df) <- colnames(h)
# ------------------------
# Heatmaps
# ------------------------
# adj.P.Val.Clade
# adj.P.Val.Salinity
# adj.P.Val.Physiology
# adj.P.Val.Salinity.Physiology
# adj.P.Val.Physiology.Clade
# adj.P.Val.Salinity.Clade
# adj.P.Val.3way

cons_salinity <- filter(new, adj.P.Val.Salinity<=0.05 & adj.P.Val.Salinity.Physiology>0.1 & adj.P.Val.3way>0.1)
h <- cons_salinity[,c(123:150)]
rownames(h)<-cons_salinity$ID
h <- h[,sample]
rownames(sample_label_df) <- colnames(h)
#--------------------------
# Set heatmap color scale
#--------------------------

my.breaks <- c(seq(-2, 0, by=0.1), seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("#00FFFF", "black"))(length(my.breaks)/2), colorRampPalette(colors = c("black", "yellow"))(length(my.breaks)/2))

annotation_colors = list(salinity = c("0ppt"="red", "15ppt"="blue"),physiology = c("M"="blue", "F"="red"),clade = c("clade1"="green", "clade2"="orange", "clade3"="purple"))

# ------------------
# Run heatmap
# ------------------

out <- pheatmap(h, cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = F,
         annotation_col = sample_label_df,
         show_rownames=F,
         cutree_rows = 2,
         color = my.colors,
         annotation_colors = annotation_colors,
         breaks = my.breaks,
         gaps_col = c(12,18)
)

# subset
cut2groups <- data.frame(sort(cutree(out$tree_row, k=2)))
colnames(cut2groups)[1] <- "group"
cons_salinity_groups <- merge(cons_salinity,cut2groups,by.x="ID",by.y="row.names")
group1up <- cons_salinity_groups[cons_salinity_groups$group == 1,]
group2down <- cons_salinity_groups[cons_salinity_groups$group == 2,]
h <- group1up[,c(123:150)]
rownames(h)<-group1up$ID
h <- h[,sample]
rownames(sample_label_df) <- colnames(h)
out <- pheatmap(h, cluster_rows = TRUE,
                clustering_distance_rows = "correlation",
                cluster_cols = F,
                annotation_col = sample_label_df,
                show_rownames=F,
                color = my.colors,
                annotation_colors = annotation_colors,
                breaks = my.breaks,
                gaps_col = c(12,18)
)
# could choose based on only pos across all 15 in h
# could choose only group1, because visually group looks up (even though not all are)
# pick a few first criteria, see how many don't match second
up <- cons_salinity[,c(123:150)]
rownames(up)<-cons_salinity$ID
# only 15
up <- up[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
up <- up[apply(up > 0, 1, all), ]
# only 14
upID <- rownames(up)
new_up <- new[new$ID %in% upID,]
new_up <- merge(new_up,cut2groups,by.x="ID",by.y="row.names")

down <- cons_salinity[,c(123:150)]
rownames(down)<-cons_salinity$ID
down <- down[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
down <- down[apply(down < 0, 1, all), ]
# only 4
downID <- rownames(down)
new_down <- new[new$ID %in% downID,]
new_down <- merge(new_down,cut2groups,by.x="ID",by.y="row.names")

# these go for GO analysis
length(group1up$ID)
length(group2down$ID)


