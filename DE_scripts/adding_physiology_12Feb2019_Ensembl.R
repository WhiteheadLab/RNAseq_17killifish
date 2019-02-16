#------------------------

# Install function for packages

#------------------------

packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
bioconductors <- function(x){
  x<- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)){
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkgs=x)
    require(x, character.only = TRUE)
  }
}
# for QC
packages("RColorBrewer")
packages("pheatmap")
packages("vsn")
packages(pheatmap)
# for DESeq

packages(lattice)
packages(RColorBrewer)
bioconductors(biomaRt)
packages(dplyr)
packages(tidyr)
packages(ggplot2)
packages(gplots)
bioconductors(DESeq2)
bioconductors(tximport)

#------------------------
library(gridExtra)
# for QA
library("RColorBrewer")
library("pheatmap")
library("vsn")
# for DESeq
library(DESeq2)
library(RColorBrewer)
library(tximport)
library(lattice)
library(gplots)
library("biomaRt")
library(dplyr)
library(ggplot2)
#library('gtable')
#library('grid')
#library('magrittr')
library(tidyr)
#------------------------


#------------------------

# Import custom scripts and data

#------------------------

if(!file.exists('plotPCAWithSampleNames.R')){
  download.file('https://gist.githubusercontent.com/ljcohen/d6cf3367efae60d7fa1ea383fd7a1296/raw/f11030ced8c7953be6fada0325b92c20d369e0a7/plotPCAWithSampleNames.R', 'plotPCAWithSampleNames.R')
}
source('plotPCAWithSampleNames.R')
if(!file.exists('overLapper_original.R')){
  download.file("https://gist.githubusercontent.com/ljcohen/b7e5fa93b8b77c33d26e4c44cadb5bb7/raw/3a2f84be2e937ef721cafaf308ff6456168f30d9/overLapper_original.R",'overLapper_original.R')
}
source('overLapper_original.R')
# This is just the counts with Experimental Design Info in the last 5 rows
if(!file.exists('../../Ensembl_species_counts_designfactors.csv')){
  download.file("https://osf.io/7vp38/download",'../../Ensembl_species_counts_designfactors.csv')
}
counts_design <- read.csv("../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)


#------------------------

# Format design and counts

#------------------------

design <- counts_design[counts_design$Ensembl == 'Empty',]
#design$type <- c("species","native_salinity","clade","group","condition")
drops <- c("X","Ensembl")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
dim(design)
dim(counts)


#------------------------

# Filter, all 128 samples have a count of at least 0.01

#------------------------
filter <- rownames(counts[rowSums(counts >= 0.01) >= 64,])
#filter <- rownames(counts[rowSums(counts >= 0.01) >= 100,])
#filter <- rownames(counts[rowSums(counts >= 0.01) >= 128,])
filtered_counts <- counts[filter,]
dim(filtered_counts)

#------------------------

# goi = Gene of Interest
# Check whether genes of interest are present after filtering
# what is the mininum amount of filtering we can get away without getting rid of these genes?

#------------------------

filtered_genes <- rownames(filtered_counts)

# ---------------------------
# Andrew Whitehead's genes of interest:
# ---------------------------

# Funhe2EKm029929
# zymogen granule membrane protein 16
goi <- filtered_genes[filtered_genes == "ENSFHEP00000007220.1"]
goi
# zymogen granule membrane protein 16
# Funhe2EKm029931
goi <- filtered_genes[filtered_genes == "ENSFHEP00000025841"]
goi
# solute carrier family 12 member 3-like (removed) 
# Funhe2EKm006896
goi <- filtered_genes[filtered_genes == "ENSFHEP00000009214"]
goi
# chloride channel, voltage-sensitive 2 (clcn2), transcript variant X2 
# Funhe2EKm024148
#goi <- res$row[res$row == "XP_012718665.1"]
goi <- filtered_genes[filtered_genes == "ENSFHEP00000019510"]
goi
# ATP-sensitive inward rectifier potassium channel 1 
# Funhe2EKm001965
goi <- filtered_genes[filtered_genes == "ENSFHEP00000015383"]
goi
# inward rectifier potassium channel 2
#Funhe2EKm023780
goi <- filtered_genes[filtered_genes == "ENSFHEP00000009753"]
goi
# --------------------------------
# other salinity genes of interest
# --------------------------------
# aquaporin-3
goi <- filtered_genes[filtered_genes == "ENSFHEP00000006725"]
goi
# cftr
# Funhe2EKm024501
goi <- filtered_genes[filtered_genes == "ENSFHEP00000008393"]
goi
# polyamine-modulated factor 1-like
# Funhe2EKm031049
goi <- filtered_genes[filtered_genes == "ENSFHEP00000013324"]
goi
# sodium/potassium/calcium exchanger 5 isoform X2
goi <- filtered_genes[filtered_genes == "ENSFHEP00000001609"]
goi
# polyamine-modulated factor 1-like
# ENSFHEP00000013324
# Funhe2EKm031049
goi <- filtered_genes[filtered_genes == "ENSFHEP00000013324"]
goi
# sodium/potassium/calcium exchanger 2
# ENSFHEP00000034177
# Funhe2EKm025497
goi <- filtered_genes[filtered_genes == "ENSFHEP00000034177"]
goi
# septin-2B isoform X2
# ENSFHEP00000015765
goi <- filtered_genes[filtered_genes == "ENSFHEP00000015765"]
goi
# CLOCK-interacting pacemaker-like
# ENSFHEP00000017303
# Funhe2EKm026846
goi <- filtered_genes[filtered_genes == "ENSFHEP00000017303"]
goi
# vasopressin V2 receptor-like
# Funhe2EKm026721
goi <- filtered_genes[filtered_genes == "ENSFHEP00000000036"]
goi
# sodium/potassium-transporting ATPase subunit beta-1-interacting protein 1
# ENSFHEP00000031108
# Funhe2EKm001174
goi <- filtered_genes[filtered_genes == "ENSFHEP00000031108"]
goi
# septin-2
# Funhe2EKm012182
goi <- filtered_genes[filtered_genes == "ENSFHEP00000016853"]
goi
# otopetrin-2
# Funhe2EKm035427
goi <- filtered_genes[filtered_genes == "ENSFHEP00000026411"]
goi
# claudin 8
# Funhe2EKm037718
goi <- filtered_genes[filtered_genes == "ENSFHEP00000006282"]
goi
# claudin 4
# ENSFHEP00000003908
# Funhe2EKm007149
goi <- filtered_genes[filtered_genes == "ENSFHEP00000003908"]
goi
# If these genes are present, then it is okay to proceed.

#------------------------

# Make Design categories

#------------------------

# design cateogories (full)
species<-as.character(unlist(design[1,]))
physiology<-as.character(unlist(design[2,]))
clade<-as.character(unlist(design[3,]))
np_cl<-as.character(unlist(design[4,]))
condition<-as.character(unlist(design[5,]))
species_physiology<-as.vector(paste(np_cl, species, sep="_"))
species_condition<-as.vector(paste(species,condition,sep="_"))
species_group<-as.vector(paste(np_cl,species,sep="_"))
species_clade<-as.vector(paste(clade,species,sep="_"))

#------------------------

# Make Design

#------------------------

cols<-colnames(counts)
ExpDesign <- data.frame(row.names=cols, physiology=physiology,condition=condition, species = species, clade=clade,species_condition=species_condition,species_clade = species_clade)
ExpDesign


m1 <- model.matrix(~condition + species_clade + species_clade:condition,ExpDesign)
colnames(m1)
all.zero <- apply(m1, 2, function(x) all(x==0))
idx <- which(all.zero)
m1 <- m1[,-idx]

all(rownames(ExpDesign) == colnames(filtered_counts))
counts_round<- round(data.matrix(filtered_counts),digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = m1)
# This takes ~50 min to run on Jetstream instance s1.xlarge (CPU: 24, Mem: 60 GB, Disk: 240 GB, Disk: 240 GB root)
# > 1 hr on Macbook Pro, 16 GB RAM

dds <- DESeq(dds, full = m1, betaPrior=FALSE)

ddsClean <- dds[which(mcols(dds)$betaConv),]
dds<-ddsClean
save(dds,file="../../dds_interactionConditionSpecies.RData")
#==========================================

# QA of DESeq

#==========================================
resultsNames(dds)
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
#ntd <- normTransform(dds)

#meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
plotDispEsts(dds)

plotPCA(vsd, intgroup=c("species_clade"))
plotPCAWithSampleNames(vsd,intgroup=c("species_clade"),ntop=40000)

# ---------------------------
# make a matrix of counts
# ---------------------------

counts_table = counts(dds, normalized=TRUE)
write.csv(counts_table,"../../killifish_counts_filtered_64samples0.01.csv")


# ---------------------------

# Get Results

# ---------------------------
resultsNames(dds)

# what are the genes that are different in BW treatment compared to FW in each species?

# Clade1_F_catanatus (reference level)
res_Fcatanatus_BW_v_FW <- results(dds, tidy=TRUE, name="condition15_ppt") %>% arrange(padj) %>% tbl_df()

# all others
res_Fdiaphanus_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_diaphanus"))) %>% arrange(padj) %>% tbl_df()
res_Fgrandis_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_grandis"))) %>% arrange(padj) %>% tbl_df()
res_FheteroclitusMDPL_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_heteroclitusMDPL"))) %>% arrange(padj) %>% tbl_df() 
res_FheteroclitusMDPP_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_heteroclitusMDPP"))) %>% arrange(padj) %>% tbl_df() 
res_Frathbuni_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_rathbuni"))) %>% arrange(padj) %>% tbl_df()
res_Fsimilis_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade1_F_similis"))) %>% arrange(padj) %>% tbl_df()         
res_Fparvapinis_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade2_F_parvapinis"))) %>% arrange(padj) %>% tbl_df()     
res_Lgoodei_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade2_L_goodei"))) %>% arrange(padj) %>% tbl_df()        
res_Lparva_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade2_L_parva"))) %>% arrange(padj) %>% tbl_df()           
res_Axenica_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_A_xenica"))) %>% arrange(padj) %>% tbl_df()       
res_Fchrysotus_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_F_chrysotus"))) %>% arrange(padj) %>% tbl_df()             
res_Fnotatus_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_F_notatus"))) %>% arrange(padj) %>% tbl_df()          
res_Folivaceous_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_F_olivaceous"  ))) %>% arrange(padj) %>% tbl_df()   
res_Fsciadicus_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_F_sciadicus"))) %>% arrange(padj) %>% tbl_df()   
res_Fzebrinus_BW_v_FW <- results(dds, tidy=TRUE, contrast=list(c("condition15_ppt","condition15_ppt.species_cladeClade3_F_zebrinus"))) %>% arrange(padj) %>% tbl_df()           

# write all to files
write.csv(res_Fcatanatus_BW_v_FW,"../../res_Fcatanatus_BW_v_FW.csv")
write.csv(res_Fdiaphanus_BW_v_FW,"../../res_Fdiaphanus_BW_v_FW.csv")
write.csv(res_Fgrandis_BW_v_FW,"../../res_Fgrandis_BW_v_FW.csv")
write.csv(res_FheteroclitusMDPL_BW_v_FW,"../../res_FheteroclitusMDPL_BW_v_FW.csv")
write.csv(res_FheteroclitusMDPP_BW_v_FW,"../../res_FheteroclitusMDPP_BW_v_FW.csv")
write.csv(res_Frathbuni_BW_v_FW,"../../res_Frathbuni_BW_v_FW.csv")
write.csv(res_Fsimilis_BW_v_FW,"../../res_Fsimilis_BW_v_FW.csv")
write.csv(res_Fparvapinis_BW_v_FW,"../../res_Fparvapinis_BW_v_FW.csv")
write.csv(res_Lgoodei_BW_v_FW,"../../res_Lgoodei_BW_v_FW.csv")
write.csv(res_Lparva_BW_v_FW,"../../res_Lparva_BW_v_FW.csv")
write.csv(res_Axenica_BW_v_FW,"../../res_Axenica_BW_v_FW.csv")
write.csv(res_Fchrysotus_BW_v_FW,"../../res_Fchrysotus_BW_v_FW.csv")
write.csv(res_Fnotatus_BW_v_FW,"../../res_Fnotatus_BW_v_FW.csv")
write.csv(res_Folivaceous_BW_v_FW,"../../res_Folivaceous_BW_v_FW.csv")
write.csv(res_Fsciadicus_BW_v_FW,"../../res_Fsciadicus_BW_v_FW.csv")
write.csv(res_Fzebrinus_BW_v_FW,"../../res_Fzebrinus_BW_v_FW.csv")
       
# what are the genes that are different in TR treatment compared to FW in each species?

# Clade1_F_catanatus (reference level)
res_Fcatanatus_TR_v_FW <- results(dds, tidy=TRUE, name="conditiontransfer") %>% arrange(padj) %>% tbl_df()

# all others
res_Fdiaphanus_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_diaphanus"))) %>% arrange(padj) %>% tbl_df()
res_Fgrandis_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_grandis"))) %>% arrange(padj) %>% tbl_df()
res_FheteroclitusMDPL_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_heteroclitusMDPL"))) %>% arrange(padj) %>% tbl_df() 
res_FheteroclitusMDPP_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_heteroclitusMDPP"))) %>% arrange(padj) %>% tbl_df() 
res_Frathbuni_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_rathbuni"))) %>% arrange(padj) %>% tbl_df()
res_Fsimilis_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade1_F_similis"))) %>% arrange(padj) %>% tbl_df()         
res_Fparvapinis_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade2_F_parvapinis"))) %>% arrange(padj) %>% tbl_df()     
res_Lgoodei_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade2_L_goodei"))) %>% arrange(padj) %>% tbl_df()        
res_Lparva_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade2_L_parva"))) %>% arrange(padj) %>% tbl_df()           
res_Axenica_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_A_xenica"))) %>% arrange(padj) %>% tbl_df()       
res_Fchrysotus_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_F_chrysotus"))) %>% arrange(padj) %>% tbl_df()             
res_Fnotatus_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_F_notatus"))) %>% arrange(padj) %>% tbl_df()          
res_Folivaceous_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_F_olivaceous"  ))) %>% arrange(padj) %>% tbl_df()   
res_Fsciadicus_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_F_sciadicus"))) %>% arrange(padj) %>% tbl_df()   
# not enough samples for this contrast in F. zebrinus
#res_Fzebrinus_TR_v_FW <- results(dds, tidy=TRUE, contrast=list(c("conditiontransfer","conditiontransfer.species_cladeClade3_F_zebrinus"))) %>% arrange(padj) %>% tbl_df()           

# write all to files
write.csv(res_Fcatanatus_TR_v_FW,"../../res_Fcatanatus_TR_v_FW.csv")
write.csv(res_Fdiaphanus_TR_v_FW,"../../res_Fdiaphanus_TR_v_FW.csv")
write.csv(res_Fgrandis_TR_v_FW,"../../res_Fgrandis_TR_v_FW.csv")
write.csv(res_FheteroclitusMDPL_TR_v_FW,"../../res_FheteroclitusMDPL_TR_v_FW.csv")
write.csv(res_FheteroclitusMDPP_TR_v_FW,"../../res_FheteroclitusMDPP_TR_v_FW.csv")
write.csv(res_Frathbuni_TR_v_FW,"../../res_Frathbuni_TR_v_FW.csv")
write.csv(res_Fsimilis_TR_v_FW,"../../res_Fsimilis_TR_v_FW.csv")
write.csv(res_Fparvapinis_TR_v_FW,"../../res_Fparvapinis_TR_v_FW.csv")
write.csv(res_Lgoodei_TR_v_FW,"../../res_Lgoodei_TR_v_FW.csv")
write.csv(res_Lparva_TR_v_FW,"../../res_Lparva_TR_v_FW.csv")
write.csv(res_Axenica_TR_v_FW,"../../res_Axenica_TR_v_FW.csv")
write.csv(res_Fchrysotus_TR_v_FW,"../../res_Fchrysotus_TR_v_FW.csv")
write.csv(res_Fnotatus_TR_v_FW,"../../res_Fnotatus_TR_v_FW.csv")
write.csv(res_Folivaceous_TR_v_FW,"../../res_Folivaceous_TR_v_FW.csv")
write.csv(res_Fsciadicus_TR_v_FW,"../../res_Fsciadicus_TR_v_FW.csv")

# ============================================
# biomart annotation!
# https://uswest.ensembl.org/Fundulus_heteroclitus/Info/Index
# ============================================

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
ensembl_proteinID = rownames(counts_table)
length(ensembl_proteinID)
#query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','go_id','description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
head(query)
dim(query)
length(unique(query$ensembl_peptide_id))
query <- query[!duplicated(query[,c(1)]),]
dim(query)
# ============================================

# Make stats tables
# one for each species

# ============================================
#
# F. catanatus
#
# -----------------------------
# counts for F. catanatus
# -----------------------------

Fcat_cols <- c("F_catanatus_BW_1.quant", "F_catanatus_BW_2.quant", "F_catanatus_BW_3.quant", 
               "F_catanatus_FW_1.quant", "F_catanatus_FW_2.quant",
               "F_catanatus_transfer_1.quant","F_catanatus_transfer_2.quant")
Fcat_counts_table <- as.data.frame(counts_table[,Fcat_cols])
dim(Fcat_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fcatanatus_TR_v_FW)[names(res_Fcatanatus_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fcatanatus_BW_v_FW)[names(res_Fcatanatus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
Fcat_counts_table_stats <- merge(as.data.frame(res_Fcatanatus_TR_v_FW),Fcat_counts_table,by=0)
Fcat_counts_table_stats <- merge(as.data.frame(res_Fcatanatus_BW_v_FW),Fcat_counts_table_stats,by='row')
head(Fcat_counts_table_stats)
dim(Fcat_counts_table_stats)
Fcat_counts_table_stats <- Fcat_counts_table_stats[,-Fcat_counts_table_stats$Row.names]
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Fcat_counts_table_ann <- merge(query,Fcat_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fcat_counts_table_ann) <- Fcat_counts_table_ann$ensembl_peptide_id
dim(Fcat_counts_table_ann)
# [1] 18356    20
Fcat_counts_table_ann <- Fcat_counts_table_ann[ , -which(names(Fcat_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. diaphanus
# ==============================
# -----------------------------
# counts for F. diaphanus
# -----------------------------

Fdia_cols <- c("F_diaphanus_BW_1.quant", "F_diaphanus_BW_2.quant", 
               "F_diaphanus_FW_2.quant", "F_diaphanus_FW_3.quant",
               "F_diaphanus_transfer_1.quant","F_diaphanus_transfer_2.quant")
Fdia_counts_table <- as.data.frame(counts_table[,Fdia_cols])
dim(Fdia_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fdiaphanus_TR_v_FW)[names(res_Fdiaphanus_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fdiaphanus_BW_v_FW)[names(res_Fdiaphanus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Fdiaphanus_TR_v_FW <- as.data.frame(res_Fdiaphanus_TR_v_FW)
rownames(res_Fdiaphanus_TR_v_FW) <- res_Fdiaphanus_TR_v_FW$row
Fdia_counts_table_stats <- merge(as.data.frame(res_Fdiaphanus_TR_v_FW),Fdia_counts_table,by=0)
Fdia_counts_table_stats <- merge(as.data.frame(res_Fdiaphanus_BW_v_FW),Fdia_counts_table_stats,by='row')
head(Fdia_counts_table_stats)
dim(Fdia_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Fdia_counts_table_ann <- merge(query,Fdia_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fdia_counts_table_ann) <- Fdia_counts_table_ann$ensembl_peptide_id
dim(Fdia_counts_table_ann)
# [1] 18356    20
Fdia_counts_table_ann <- Fdia_counts_table_ann[ , -which(names(Fdia_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. grandis
# ==============================
# -----------------------------
# counts for F. grandis
# -----------------------------

Fgran_cols <- c("F_grandis_BW_1.quant", "F_grandis_BW_2.quant", "F_grandis_BW_3.quant",
               "F_grandis_FW_1.quant", "F_grandis_FW_2.quant", "F_grandis_FW_3.quant",
               "F_grandis_transfer_1.quant", "F_grandis_transfer_2.quant", "F_grandis_transfer_3.quant")
Fgran_counts_table <- as.data.frame(counts_table[,Fgran_cols])
dim(Fgran_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fgrandis_TR_v_FW)[names(res_Fgrandis_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fgrandis_BW_v_FW)[names(res_Fgrandis_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Fgrandis_TR_v_FW <- as.data.frame(res_Fgrandis_TR_v_FW)
rownames(res_Fgrandis_TR_v_FW) <- res_Fgrandis_TR_v_FW$row
Fgran_counts_table_stats <- merge(as.data.frame(res_Fgrandis_TR_v_FW),Fgran_counts_table,by=0)
Fgran_counts_table_stats <- merge(as.data.frame(res_Fgrandis_BW_v_FW),Fgran_counts_table_stats,by='row')
head(Fgran_counts_table_stats)
dim(Fgran_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Fgran_counts_table_ann <- merge(query,Fgran_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fgran_counts_table_ann) <- Fgran_counts_table_ann$ensembl_peptide_id
dim(Fgran_counts_table_ann)
# [1] 18356    20
Fgran_counts_table_ann <- Fgran_counts_table_ann[ , -which(names(Fgran_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. heteroclitus (MDPL)
# ==============================
# -----------------------------
# counts for F. heteroclitus (MDPL)
# -----------------------------

FhetMDPL_cols <- c("F_heteroclitusMDPL_BW_1.quant" , "F_heteroclitusMDPL_BW_2.quant" , "F_heteroclitusMDPL_BW_3.quant" ,
                "F_heteroclitusMDPL_FW_1.quant", "F_heteroclitusMDPL_FW_2.quant", "F_heteroclitusMDPL_FW_3.quant" ,
                "F_heteroclitusMDPL_transfer_1.quant", "F_heteroclitusMDPL_transfer_2.quant", "F_heteroclitusMDPL_transfer_3.quant")
FhetMDPL_counts_table <- as.data.frame(counts_table[,FhetMDPL_cols])
dim(FhetMDPL_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_FheteroclitusMDPL_TR_v_FW)[names(res_FheteroclitusMDPL_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPL_BW_v_FW)[names(res_FheteroclitusMDPL_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_FheteroclitusMDPL_TR_v_FW <- as.data.frame(res_FheteroclitusMDPL_TR_v_FW)
rownames(res_FheteroclitusMDPL_TR_v_FW) <- res_FheteroclitusMDPL_TR_v_FW$row
FhetMDPL_counts_table_stats <- merge(as.data.frame(res_FheteroclitusMDPL_TR_v_FW),FhetMDPL_counts_table,by=0)
FhetMDPL_counts_table_stats <- merge(as.data.frame(res_FheteroclitusMDPL_BW_v_FW),FhetMDPL_counts_table_stats,by='row')
head(FhetMDPL_counts_table_stats)
dim(FhetMDPL_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
FhetMDPL_counts_table_ann <- merge(query,FhetMDPL_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(FhetMDPL_counts_table_ann) <- FhetMDPL_counts_table_ann$ensembl_peptide_id
dim(FhetMDPL_counts_table_ann)
# [1] 18356    20
FhetMDPL_counts_table_ann <- FhetMDPL_counts_table_ann[ , -which(names(FhetMDPL_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. heteroclitus (MDPP)
# ==============================
# -----------------------------
# counts for F. heteroclitus (MDPP)
# -----------------------------

FhetMDPP_cols <- c("F_heteroclitusMDPP_BW_1.quant" , "F_heteroclitusMDPP_BW_2.quant" , "F_heteroclitusMDPP_BW_3.quant" ,
                   "F_heteroclitusMDPP_FW_1.quant", "F_heteroclitusMDPP_FW_2.quant", "F_heteroclitusMDPP_FW_3.quant" ,
                   "F_heteroclitusMDPP_transfer_1.quant", "F_heteroclitusMDPP_transfer_2.quant", "F_heteroclitusMDPP_transfer_3.quant")
FhetMDPP_counts_table <- as.data.frame(counts_table[,FhetMDPP_cols])
dim(FhetMDPP_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_FheteroclitusMDPP_TR_v_FW)[names(res_FheteroclitusMDPP_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_FheteroclitusMDPP_BW_v_FW)[names(res_FheteroclitusMDPP_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_FheteroclitusMDPP_TR_v_FW <- as.data.frame(res_FheteroclitusMDPP_TR_v_FW)
rownames(res_FheteroclitusMDPP_TR_v_FW) <- res_FheteroclitusMDPP_TR_v_FW$row
FhetMDPP_counts_table_stats <- merge(as.data.frame(res_FheteroclitusMDPP_TR_v_FW),FhetMDPP_counts_table,by=0)
FhetMDPP_counts_table_stats <- merge(as.data.frame(res_FheteroclitusMDPP_BW_v_FW),FhetMDPP_counts_table_stats,by='row')
head(FhetMDPP_counts_table_stats)
dim(FhetMDPP_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
FhetMDPP_counts_table_ann <- merge(query,FhetMDPP_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(FhetMDPP_counts_table_ann) <- FhetMDPP_counts_table_ann$ensembl_peptide_id
dim(FhetMDPP_counts_table_ann)
# [1] 18356    20
FhetMDPP_counts_table_ann <- FhetMDPP_counts_table_ann[ , -which(names(FhetMDPP_counts_table_ann) %in% c("Row.names"))]


# ==============================
# F. rathbuni
# ==============================
# -----------------------------
# counts for F. heteroclitus (MDPP)
# -----------------------------

Frath_cols <- c("F_rathbuni_BW_1.quant","F_rathbuni_BW_2.quant","F_rathbuni_BW_3.quant",
                   "F_rathbuni_FW_1.quant","F_rathbuni_BW_2.quant","F_rathbuni_FW_3.quant",
                   "F_rathbuni_transfer_1.quant","F_rathbuni_transfer_2.quant","F_rathbuni_transfer_3.quant")
Frath_counts_table <- as.data.frame(counts_table[,Frath_cols])
dim(Frath_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Frathbuni_TR_v_FW)[names(res_Frathbuni_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Frathbuni_BW_v_FW)[names(res_Frathbuni_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Frathbuni_TR_v_FW <- as.data.frame(res_Frathbuni_TR_v_FW)
rownames(res_Frathbuni_TR_v_FW) <- res_Frathbuni_TR_v_FW$row
Frath_counts_table_stats <- merge(as.data.frame(res_Frathbuni_TR_v_FW),Frath_counts_table,by=0)
Frath_counts_table_stats <- merge(as.data.frame(res_Frathbuni_BW_v_FW),Frath_counts_table_stats,by='row')
head(Frath_counts_table_stats)
dim(Frath_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Frath_counts_table_ann <- merge(query,Frath_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Frath_counts_table_ann) <- Frath_counts_table_ann$ensembl_peptide_id
dim(FhetMDPP_counts_table_ann)
# [1] 18356    20
Frath_counts_table_ann <- Frath_counts_table_ann[ , -which(names(Frath_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. similis
# ==============================
# -----------------------------
# counts for F. simils
# -----------------------------

Fsim_cols <- c("F_similis_BW_1.quant","F_similis_BW_2.quant","F_similis_BW_3.quant",
                "F_similis_FW_1.quant","F_similis_FW_2.quant","F_similis_FW_3.quant",
                "F_similis_transfer_1.quant","F_similis_transfer_2.quant","F_similis_transfer_3.quant")
Fsim_counts_table <- as.data.frame(counts_table[,Fsim_cols])
dim(Fsim_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fsimilis_TR_v_FW)[names(res_Fsimilis_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fsimilis_BW_v_FW)[names(res_Fsimilis_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Fsimilis_TR_v_FW <- as.data.frame(res_Fsimilis_TR_v_FW)
rownames(res_Fsimilis_TR_v_FW) <- res_Fsimilis_TR_v_FW$row
Fsim_counts_table_stats <- merge(as.data.frame(res_Fsimilis_TR_v_FW),Fsim_counts_table,by=0)
Fsim_counts_table_stats <- merge(as.data.frame(res_Fsimilis_BW_v_FW),Fsim_counts_table_stats,by='row')
head(Fsim_counts_table_stats)
dim(Fsim_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Fsim_counts_table_ann <- merge(query,Fsim_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fsim_counts_table_ann) <- Fsim_counts_table_ann$ensembl_peptide_id
dim(Fsim_counts_table_ann)
# [1] 18356    20
Fsim_counts_table_ann <- Fsim_counts_table_ann[ , -which(names(Fsim_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. parvapinis
# ==============================
# -----------------------------
# counts for F. parvapinis
# -----------------------------

Fparv_cols <- c("F_parvapinis_BW_1.quant","F_parvapinis_BW_2.quant","F_parvapinis_BW_3.quant",
               "F_parvapinis_FW_1.quant","F_parvapinis_FW_2.quant","F_parvapinis_FW_3.quant",
               "F_parvapinis_transfer_1.quant","F_parvapinis_transfer_2.quant")
Fparv_counts_table <- as.data.frame(counts_table[,Fparv_cols])
dim(Fparv_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fparvapinis_TR_v_FW)[names(res_Fparvapinis_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fparvapinis_BW_v_FW)[names(res_Fparvapinis_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Fparvapinis_TR_v_FW <- as.data.frame(res_Fparvapinis_TR_v_FW)
rownames(res_Fparvapinis_TR_v_FW) <- res_Fparvapinis_TR_v_FW$row
Fparv_counts_table_stats <- merge(as.data.frame(res_Fparvapinis_TR_v_FW),Fparv_counts_table,by=0)
Fparv_counts_table_stats <- merge(as.data.frame(res_Fparvapinis_BW_v_FW),Fparv_counts_table_stats,by='row')
head(Fparv_counts_table_stats)
dim(Fparv_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Fparv_counts_table_ann <- merge(query,Fparv_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fparv_counts_table_ann) <- Fparv_counts_table_ann$ensembl_peptide_id
dim(Fparv_counts_table_ann)
# [1] 18356    20
Fparv_counts_table_ann <- Fparv_counts_table_ann[ , -which(names(Fparv_counts_table_ann) %in% c("Row.names"))]

# ==============================
# L. goodei
# ==============================
# -----------------------------
# counts for L. goodei
# -----------------------------

Lgod_cols <- c("L_goodei_BW_1.quant","L_goodei_BW_2.quant","L_goodei_BW_3.quant",
                "L_goodei_FW_1.quant","L_goodei_FW_2.quant","L_goodei_FW_3.quant",
                "L_goodei_transfer_1.quant","L_goodei_transfer_2.quant","L_goodei_transfer_3.quant")
Lgod_counts_table <- as.data.frame(counts_table[,Lgod_cols])
dim(Lgod_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Lgoodei_TR_v_FW)[names(res_Lgoodei_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Lgoodei_BW_v_FW)[names(res_Lgoodei_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and TR_v_FW stats
# -----------------------------
res_Lgoodei_TR_v_FW <- as.data.frame(res_Lgoodei_TR_v_FW)
rownames(res_Lgoodei_TR_v_FW) <- res_Lgoodei_TR_v_FW$row
Lgod_counts_table_stats <- merge(as.data.frame(res_Lgoodei_TR_v_FW),Lgod_counts_table,by=0)
Lgod_counts_table_stats <- merge(as.data.frame(res_Lgoodei_BW_v_FW),Lgod_counts_table_stats,by='row')
head(Lgod_counts_table_stats)
dim(Lgod_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Lgod_counts_table_ann <- merge(query,Lgod_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Lgod_counts_table_ann) <- Lgod_counts_table_ann$ensembl_peptide_id
dim(Lgod_counts_table_ann)
# [1] 18356    20
Lgod_counts_table_ann <- Lgod_counts_table_ann[ , -which(names(Lgod_counts_table_ann) %in% c("Row.names"))]

# ==============================
# L. parva
# ==============================
# -----------------------------
# counts for L. parva
# -----------------------------

Lpar_cols <- c("L_parva_BW_1.quant","L_parva_BW_2.quant","L_parva_BW_3.quant",
               "L_parva_FW_1.quant","L_parva_FW_2.quant","L_parva_FW_3.quant",
               "L_parva_transfer_1.quant","L_parva_transfer_2.quant","L_parva_transfer_3.quant")
Lpar_counts_table <- as.data.frame(counts_table[,Lpar_cols])
dim(Lpar_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Lparva_TR_v_FW)[names(res_Lparva_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Lparva_BW_v_FW)[names(res_Lparva_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Lparva_TR_v_FW <- as.data.frame(res_Lparva_TR_v_FW)
rownames(res_Lparva_TR_v_FW) <- res_Lparva_TR_v_FW$row
Lpar_counts_table_stats <- merge(as.data.frame(res_Lparva_TR_v_FW),Lpar_counts_table,by=0)
Lpar_counts_table_stats <- merge(as.data.frame(res_Lparva_BW_v_FW),Lpar_counts_table_stats,by='row')
head(Lpar_counts_table_stats)
dim(Lpar_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Lpar_counts_table_ann <- merge(query,Lpar_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Lpar_counts_table_ann) <- Lpar_counts_table_ann$ensembl_peptide_id
dim(Lpar_counts_table_ann)
# [1] 18356    20
Lpar_counts_table_ann <- Lpar_counts_table_ann[ , -which(names(Lpar_counts_table_ann) %in% c("Row.names"))]

# ==============================
# A. xenica
# ==============================
# -----------------------------
# counts for A. xenica
# -----------------------------

Axen_cols <- c("A_xenica_BW_1.quant","A_xenica_BW_2.quant","A_xenica_BW_3.quant",
               "A_xenica_FW_1.quant","A_xenica_FW_2.quant","A_xenica_FW_3.quant",
               "A_xenica_transfer_1.quant","A_xenica_transfer_2.quant","A_xenica_transfer_3.quant")
Axen_counts_table <- as.data.frame(counts_table[,Axen_cols])
dim(Axen_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Axenica_TR_v_FW)[names(res_Axenica_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Axenica_BW_v_FW)[names(res_Axenica_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Axenica_TR_v_FW <- as.data.frame(res_Axenica_TR_v_FW)
rownames(res_Axenica_TR_v_FW) <- res_Axenica_TR_v_FW$row
Axen_counts_table_stats <- merge(as.data.frame(res_Axenica_TR_v_FW),Axen_counts_table,by=0)
Axen_counts_table_stats <- merge(as.data.frame(res_Axenica_BW_v_FW),Axen_counts_table_stats,by='row')
head(Axen_counts_table_stats)
dim(Axen_counts_table_stats)
# [1] 18356    15
# -----------------------------
# merge annotations with stats
# -----------------------------
Axen_counts_table_ann <- merge(query,Axen_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Axen_counts_table_ann) <- Axen_counts_table_ann$ensembl_peptide_id
dim(Axen_counts_table_ann)
# [1] 18356    20
Axen_counts_table_ann <- Axen_counts_table_ann[ , -which(names(Axen_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. chrysotus
# ==============================
# -----------------------------
# counts for F. chrysotus
# -----------------------------

Fchry_cols <- c("F_chrysotus_BW_1.quant","F_chrysotus_BW_2.quant","F_chrysotus_BW_3.quant",
               "F_chrysotus_FW_1.quant","F_chrysotus_FW_2.quant","F_chrysotus_FW_3.quant",
               "F_chrysotus_transfer_1.quant","F_chrysotus_transfer_2.quant")
Fchry_counts_table <- as.data.frame(counts_table[,Fchry_cols])
dim(Fchry_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fchrysotus_TR_v_FW)[names(res_Fchrysotus_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fchrysotus_BW_v_FW)[names(res_Fchrysotus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Fchrysotus_TR_v_FW <- as.data.frame(res_Fchrysotus_TR_v_FW)
rownames(res_Fchrysotus_TR_v_FW) <- res_Fchrysotus_TR_v_FW$row
Fchry_counts_table_stats <- merge(as.data.frame(res_Fchrysotus_TR_v_FW),Fchry_counts_table,by=0)
Fchry_counts_table_stats <- merge(as.data.frame(res_Fchrysotus_BW_v_FW),Fchry_counts_table_stats,by='row')
head(Fchry_counts_table_stats)
dim(Fchry_counts_table_stats)
# [1] 18356    15

# -----------------------------
# merge annotations with stats
# -----------------------------
Fchry_counts_table_ann <- merge(query,Fchry_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fchry_counts_table_ann) <- Fchry_counts_table_ann$ensembl_peptide_id
dim(Fchry_counts_table_ann)
# [1] 18356    20
Fchry_counts_table_ann <- Fchry_counts_table_ann[ , -which(names(Fchry_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. notatus
# ==============================
# -----------------------------
# counts for F. notatus
# -----------------------------

Fnota_cols <- c("F_notatus_BW_1.quant","F_notatus_BW_2.quant","F_notatus_BW_3.quant",
                "F_notatus_FW_1.quant","F_notatus_FW_2.quant","F_notatus_FW_3.quant",
                "F_notatus_transfer_1.quant","F_notatus_transfer_2.quant","F_notatus_transfer_3.quant")
Fnota_counts_table <- as.data.frame(counts_table[,Fnota_cols])
dim(Fnota_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fnotatus_TR_v_FW)[names(res_Fnotatus_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fnotatus_BW_v_FW)[names(res_Fnotatus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Fnotatus_TR_v_FW <- as.data.frame(res_Fnotatus_TR_v_FW)
rownames(res_Fnotatus_TR_v_FW) <- res_Fnotatus_TR_v_FW$row
Fnota_counts_table_stats <- merge(as.data.frame(res_Fnotatus_TR_v_FW),Fnota_counts_table,by=0)
Fnota_counts_table_stats <- merge(as.data.frame(res_Fnotatus_BW_v_FW),Fnota_counts_table_stats,by='row')
head(Fnota_counts_table_stats)
dim(Fnota_counts_table_stats)
# [1] 18356    15

# -----------------------------
# merge annotations with stats
# -----------------------------
Fnota_counts_table_ann <- merge(query,Fnota_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fnota_counts_table_ann) <- Fnota_counts_table_ann$ensembl_peptide_id
dim(Fnota_counts_table_ann)
# [1] 18356    20
Fnota_counts_table_ann <- Fnota_counts_table_ann[ , -which(names(Fnota_counts_table_ann) %in% c("Row.names"))]


# ==============================
# F. olivaceus
# ==============================
# -----------------------------
# counts for F. olivaceus
# -----------------------------

Foli_cols <- c("F_olivaceous_BW_1.quant","F_olivaceous_BW_2.quant","F_olivaceous_BW_3.quant",
                "F_olivaceous_FW_1.quant","F_olivaceous_FW_2.quant","F_olivaceous_FW_3.quant",
                "F_olivaceous_transfer_1.quant","F_olivaceous_transfer_2.quant")
Foli_counts_table <- as.data.frame(counts_table[,Foli_cols])
dim(Foli_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Folivaceous_TR_v_FW)[names(res_Folivaceous_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Folivaceous_BW_v_FW)[names(res_Folivaceous_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Folivaceous_TR_v_FW <- as.data.frame(res_Folivaceous_TR_v_FW)
rownames(res_Folivaceous_TR_v_FW) <- res_Folivaceous_TR_v_FW$row
Foli_counts_table_stats <- merge(as.data.frame(res_Folivaceous_TR_v_FW),Foli_counts_table,by=0)
Foli_counts_table_stats <- merge(as.data.frame(res_Folivaceous_BW_v_FW),Foli_counts_table_stats,by='row')
head(Foli_counts_table_stats)
dim(Foli_counts_table_stats)
# [1] 18356    15

# -----------------------------
# merge annotations with stats
# -----------------------------
Foli_counts_table_ann <- merge(query,Foli_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Foli_counts_table_ann) <- Foli_counts_table_ann$ensembl_peptide_id
dim(Foli_counts_table_ann)
# [1] 18356    20
Foli_counts_table_ann <- Foli_counts_table_ann[ , -which(names(Foli_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. sciadicus
# ==============================
# -----------------------------
# counts for F. sciadicus
# -----------------------------

Fsci_cols <- c("F_sciadicus_BW_1.quant",
               "F_sciadicus_FW_1.quant","F_sciadicus_FW_2.quant",
               "F_sciadicus_transfer_1.quant")
Fsci_counts_table <- as.data.frame(counts_table[,Fsci_cols])
dim(Fsci_counts_table)

# -----------------------------
# column names for stats from TR_v_FW specific contrast
# -----------------------------
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'padj'] <- 'padj-TR-v-0.2ppt'
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-TR-v-0.2ppt'
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'lfcSE'] <- 'lfcSE-TR-v-0.2ppt'
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'stat'] <- 'stat-TR-v-0.2ppt'
names(res_Fsciadicus_TR_v_FW)[names(res_Fsciadicus_TR_v_FW) == 'pvalue'] <- 'pvalue-TR-v-0.2ppt'
# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fsciadicus_BW_v_FW)[names(res_Fsciadicus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Fsciadicus_TR_v_FW <- as.data.frame(res_Fsciadicus_TR_v_FW)
rownames(res_Fsciadicus_TR_v_FW) <- res_Fsciadicus_TR_v_FW$row
Fsci_counts_table_stats <- merge(as.data.frame(res_Fsciadicus_TR_v_FW),Fsci_counts_table,by=0)
Fsci_counts_table_stats <- merge(as.data.frame(res_Fsciadicus_BW_v_FW),Fsci_counts_table_stats,by='row')
head(Fsci_counts_table_stats)
dim(Fsci_counts_table_stats)
# [1] 18356    15

# -----------------------------
# merge annotations with stats
# -----------------------------
Fsci_counts_table_ann <- merge(query,Fsci_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fsci_counts_table_ann) <- Fsci_counts_table_ann$ensembl_peptide_id
dim(Fsci_counts_table_ann)
# [1] 18356    20
Fsci_counts_table_ann <- Fsci_counts_table_ann[ , -which(names(Fsci_counts_table_ann) %in% c("Row.names"))]

# ==============================
# F. zebrinus
# ==============================
# -----------------------------
# counts for F. zebrinus
# -----------------------------

Fzeb_cols <- c("F_zebrinus_BW_1.quant","F_zebrinus_BW_2.quant",
               "F_zebrinus_FW_1.quant","F_zebrinus_FW_2.quant")
Fzeb_counts_table <- as.data.frame(counts_table[,Fzeb_cols])
dim(Fzeb_counts_table)

# -----------------------------
# column names for stats from BW_v_FW specific contrast
# -----------------------------
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'padj'] <- 'padj-15ppt-v-0.2ppt'
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'baseMean'] <- 'baseMean-ALL'
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'log2FoldChange'] <- 'log2FoldChange-15ppt-v-0.2ppt'
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'lfcSE'] <- 'lfcSE-15ppt-v-0.2ppt'
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'stat'] <- 'stat-15ppt-v-0.2ppt'
names(res_Fzebrinus_BW_v_FW)[names(res_Fzebrinus_BW_v_FW) == 'pvalue'] <- 'pvalue-15ppt-v-0.2ppt'
# -----------------------------
# merge counts and stats
# -----------------------------
res_Fzebrinus_BW_v_FW <- as.data.frame(res_Fzebrinus_BW_v_FW)
rownames(res_Fzebrinus_BW_v_FW) <- res_Fzebrinus_BW_v_FW$row
Fzeb_counts_table_stats <- merge(as.data.frame(res_Fzebrinus_BW_v_FW),Fzeb_counts_table,by=0)
head(Fzeb_counts_table_stats)
dim(Fsci_counts_table_stats)
# [1] 18356    15

# -----------------------------
# merge annotations with stats
# -----------------------------
Fzeb_counts_table_ann <- merge(query,Fzeb_counts_table_stats,by.x = "ensembl_peptide_id", by.y = "row", all = TRUE)
rownames(Fzeb_counts_table_ann) <- Fzeb_counts_table_ann$ensembl_peptide_id
dim(Fzeb_counts_table_ann)
# [1] 18356    20
Fzeb_counts_table_ann <- Fzeb_counts_table_ann[ , -which(names(Fzeb_counts_table_ann) %in% c("Row.names"))]

# -----------------------------
# write csv files
# -----------------------------
write.csv(Fcat_counts_table_ann,"../../counts_stats_byspecies/Fcat_stats_annotations_counts.csv")
write.csv(Fdia_counts_table_ann,"../../counts_stats_byspecies/Fdiaphanus_stats_annotations_counts.csv")
write.csv(Fgran_counts_table_ann,"../../counts_stats_byspecies/Fgrandis_stats_annotations_counts.csv")
write.csv(FhetMDPL_counts_table_ann,"../../counts_stats_byspecies/FheteroclitusMDPL_stats_annotations_counts.csv")
write.csv(FhetMDPP_counts_table_ann,"../../counts_stats_byspecies/FheteroclitusMDPP_stats_annotations_counts.csv")
write.csv(Frath_counts_table_ann,"../../counts_stats_byspecies/Frathbuni_stats_annotations_counts.csv")
write.csv(Fsim_counts_table_ann,"../../counts_stats_byspecies/Fsimils_stats_annotations_counts.csv")
write.csv(Fparv_counts_table_ann,"../../counts_stats_byspecies/Fparvapinis_stats_annotations_counts.csv")
write.csv(Lgod_counts_table_ann,"../../counts_stats_byspecies/Lgoodei_stats_annotations_counts.csv")
write.csv(Lpar_counts_table_ann,"../../counts_stats_byspecies/Lparva_stats_annotations_counts.csv")
write.csv(Axen_counts_table_ann,"../../counts_stats_byspecies/Axenica_stats_annotations_counts.csv")
write.csv(Fchry_counts_table_ann,"../../counts_stats_byspecies/Fchrysotus_stats_annotations_counts.csv")
write.csv(Fnota_counts_table_ann,"../../counts_stats_byspecies/Fnotatus_stats_annotations_counts.csv")
write.csv(Foli_counts_table_ann,"../../counts_stats_byspecies/Folivaceus_stats_annotations_counts.csv")
write.csv(Fsci_counts_table_ann,"../../counts_stats_byspecies/Fsciadicus_stats_annotations_counts.csv")
write.csv(Fzeb_counts_table_ann,"../../counts_stats_byspecies/Fzebrinus_stats_annotations_counts.csv")
# -----------------------------
# find missing annotations
# -----------------------------

counts_table <- cbind(counts_table,ensembl_proteinID)
dim(counts_table)
length(unique(rownames(counts_table)))
colnames(counts_table)[129] <- "ensembl_peptide_id"
counts_table_ann <- merge(counts_table,query,by="ensembl_peptide_id")
counts_id <- counts_table[,129]
ann_id <- query$ensembl_peptide_id
counts_minus_ann<- counts_table[!(counts_table[,129] %in% ann_id),]
dim(counts_minus_ann)
head(counts_minus_ann)
missing_ensembl_ID <- rownames(counts_minus_ann)

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
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000007220.1"]
goi
# zymogen granule membrane protein 16
# Funhe2EKm029931
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000025841"]
goi
# solute carrier family 12 member 3-like (removed) 
# Funhe2EKm006896
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000009214"]
goi
# chloride channel, voltage-sensitive 2 (clcn2), transcript variant X2 (removed)
# Funhe2EKm024148
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000019510"]
goi
# ATP-sensitive inward rectifier potassium channel 1 
# Funhe2EKm001965
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000015383"]
goi
# inward rectifier potassium channel 2
#Funhe2EKm023780
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000009753"]

# --------------------------------
# other salinity genes of interest
# --------------------------------
# ============================================
# aquaporin-3
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000006725"]
goi
# cftr
# Funhe2EKm024501
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000008393"]
goi
# polyamine-modulated factor 1-like
# Funhe2EKm031049
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000013324"]
goi
# sodium/potassium/calcium exchanger 2
# ENSFHEP00000034177
# Funhe2EKm025497
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000034177"]
goi
# septin-2B isoform X2
# ENSFHEP00000015765
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000015765"]
goi
# CLOCK-interacting pacemaker-like
# ENSFHEP00000017303
# Funhe2EKm026846
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000017303"]
goi
# vasopressin V2 receptor-like
# Funhe2EKm026721
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000000036"]
goi
# sodium/potassium-transporting ATPase subunit beta-1-interacting protein 1
# ENSFHEP00000031108
# Funhe2EKm001174
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000031108"]
goi
# septin-2
# Funhe2EKm012182
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000016853"]
goi
# otopetrin-2
# Funhe2EKm035427
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000026411"]
goi
# claudin 8
# Funhe2EKm037718
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000006282"]
goi
# claudin 4
# ENSFHEP00000003908
# Funhe2EKm007149
goi <- res_Fdiaphanus_BW_v_FW$row[res_Fdiaphanus_BW_v_FW$row == "ENSFHEP00000003908"]
goi
all_goi<-c("ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
           "ENSFHEP00000015383","ENSFHEP00000009753","ENSFHEP00000006725","ENSFHEP00000008393",
           "ENSFHEP00000013324","ENSFHEP00000001609","ENSFHEP00000013324","ENSFHEP00000034177",
           "ENSFHEP00000015765","ENSFHEP00000017303","ENSFHEP00000000036","ENSFHEP00000031108",
           "ENSFHEP00000016853","ENSFHEP00000003908")

# ============================================
#
# Make plots with goi
# This was very helpful:
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf
# 
# ============================================

pdf("../../multi-ggplot2-catalog_salinity_14Feb2018.pdf",paper="USr",width=13.5, height=8)
for (i in all_goi){
  tcounts <- t(log2((counts(dds[i, ], normalized=TRUE, replaced=FALSE)+.5))) %>% 
    merge(colData(dds), ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(i)+1):ncol(.))
  #tcounts %>% select(Row.names, species, clade, condition, gene, expression) %>% head %>% knitr::kable()
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","transfer","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene~species,scales='fixed') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","transfer","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene~species,scales='fixed') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","transfer","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene~species,scales='fixed',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()

# ============================================
#
# Find sig genes across species
# up FC
# BW_v_FW
# Clade 1
# ============================================
# subset each stats table for sig genes
# Clade 1

Fcat_counts_table_ann_sig <- subset(Fcat_counts_table_ann,Fcat_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fcat_counts_table_ann_up <- subset(Fcat_counts_table_ann_sig,Fcat_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fcat_counts_table_ann_down <- subset(Fcat_counts_table_ann_sig,Fcat_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fcat_counts_table_ann_sig)
Fdia_counts_table_ann_sig <- subset(Fdia_counts_table_ann,Fdia_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fdia_counts_table_ann_up <- subset(Fdia_counts_table_ann_sig,Fdia_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fdia_counts_table_ann_down <- subset(Fdia_counts_table_ann_sig,Fdia_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fdia_counts_table_ann_sig)
Fgran_counts_table_ann_sig <- subset(Fgran_counts_table_ann,Fgran_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fgran_counts_table_ann_up <- subset(Fgran_counts_table_ann_sig,Fgran_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fgran_counts_table_ann_down <- subset(Fgran_counts_table_ann_sig,Fgran_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fgran_counts_table_ann_sig)
FhetMDPL_counts_table_ann_sig <- subset(FhetMDPL_counts_table_ann,FhetMDPL_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
FhetMDPL_counts_table_ann_up <- subset(FhetMDPL_counts_table_ann_sig,FhetMDPL_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
FhetMDPL_counts_table_ann_down <- subset(FhetMDPL_counts_table_ann_sig,FhetMDPL_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(FhetMDPL_counts_table_ann_sig)
FhetMDPP_counts_table_ann_sig <- subset(FhetMDPP_counts_table_ann,FhetMDPP_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
FhetMDPP_counts_table_ann_up <- subset(FhetMDPP_counts_table_ann_sig,FhetMDPP_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
FhetMDPP_counts_table_ann_down <- subset(FhetMDPP_counts_table_ann_sig,FhetMDPP_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(FhetMDPP_counts_table_ann_sig)
Frath_counts_table_ann_sig <- subset(Frath_counts_table_ann,Frath_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Frath_counts_table_ann_up <- subset(Frath_counts_table_ann_sig,Frath_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Frath_counts_table_ann_down <- subset(Frath_counts_table_ann_sig,Frath_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Frath_counts_table_ann_sig)
Fsim_counts_table_ann_sig <- subset(Fsim_counts_table_ann,Fsim_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fsim_counts_table_ann_up <- subset(Fsim_counts_table_ann_sig,Fsim_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fsim_counts_table_ann_down <- subset(Fsim_counts_table_ann_sig,Fsim_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fsim_counts_table_ann_sig)

# Clade 2
Fparv_counts_table_ann_sig <- subset(Fparv_counts_table_ann,Fparv_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fparv_counts_table_ann_up <- subset(Fparv_counts_table_ann_sig,Fparv_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fparv_counts_table_ann_down <- subset(Fparv_counts_table_ann_sig,Fparv_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fparv_counts_table_ann_sig)
Lgod_counts_table_ann_sig <- subset(Lgod_counts_table_ann,Lgod_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Lgod_counts_table_ann_up <- subset(Lgod_counts_table_ann_sig,Lgod_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Lgod_counts_table_ann_down <- subset(Lgod_counts_table_ann_sig,Lgod_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Lgod_counts_table_ann_sig)
Lpar_counts_table_ann_sig <- subset(Lpar_counts_table_ann,Lpar_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Lpar_counts_table_ann_up <- subset(Lpar_counts_table_ann_sig,Lpar_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Lpar_counts_table_ann_down <- subset(Lpar_counts_table_ann_sig,Lpar_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Lpar_counts_table_ann_sig)
# Clade 3

Axen_counts_table_ann_sig <- subset(Axen_counts_table_ann,Axen_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Axen_counts_table_ann_up <- subset(Axen_counts_table_ann_sig,Axen_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Axen_counts_table_ann_down <- subset(Axen_counts_table_ann_sig,Axen_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Axen_counts_table_ann_sig)
Fchry_counts_table_ann_sig <- subset(Fchry_counts_table_ann,Fchry_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fchry_counts_table_ann_up <- subset(Fchry_counts_table_ann_sig,Fchry_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fchry_counts_table_ann_down <- subset(Fchry_counts_table_ann_sig,Fchry_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fchry_counts_table_ann_sig)
Fnota_counts_table_ann_sig <- subset(Fnota_counts_table_ann,Fnota_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fnota_counts_table_ann_up <- subset(Fnota_counts_table_ann_sig,Fnota_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fnota_counts_table_ann_down <- subset(Fnota_counts_table_ann_sig,Fnota_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fnota_counts_table_ann_sig)
Foli_counts_table_ann_sig <- subset(Foli_counts_table_ann,Foli_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Foli_counts_table_ann_up <- subset(Foli_counts_table_ann_sig,Foli_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Foli_counts_table_ann_down <- subset(Foli_counts_table_ann_sig,Foli_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Foli_counts_table_ann_sig)
Fsci_counts_table_ann_sig <- subset(Fsci_counts_table_ann,Fsci_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fsci_counts_table_ann_up <- subset(Fsci_counts_table_ann_sig,Fsci_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fsci_counts_table_ann_down <- subset(Fsci_counts_table_ann_sig,Fsci_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fsci_counts_table_ann_sig)
Fzeb_counts_table_ann_sig <- subset(Fzeb_counts_table_ann,Fzeb_counts_table_ann$`padj-15ppt-v-0.2ppt`<0.05)
Fzeb_counts_table_ann_up <- subset(Fzeb_counts_table_ann_sig,Fzeb_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` > 0.5)
Fzeb_counts_table_ann_down <- subset(Fzeb_counts_table_ann_sig,Fzeb_counts_table_ann_sig$`log2FoldChange-15ppt-v-0.2ppt` < -0.5)
dim(Fzeb_counts_table_ann_sig)

# merge common genes up
# Clade 1
length(as.vector(Fcat_counts_table_ann_up$ensembl_peptide_id))
length(as.vector(Fdia_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(as.vector(Fcat_counts_table_ann_up$ensembl_peptide_id),as.vector(Fdia_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
length(as.vector(Fgran_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(clade1_common_genes,as.vector(Fgran_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
length(as.vector(FhetMDPL_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(clade1_common_genes,as.vector(FhetMDPL_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
length(as.vector(FhetMDPP_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(clade1_common_genes,as.vector(FhetMDPP_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
length(as.vector(Frath_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(clade1_common_genes,as.vector(Frath_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
length(as.vector(Fsim_counts_table_ann_up$ensembl_peptide_id))
clade1_common_genes <- union(clade1_common_genes,as.vector(Fsim_counts_table_ann_up$ensembl_peptide_id))
length(clade1_common_genes)
# Clade 2
length(as.vector(Fparv_counts_table_ann_up$ensembl_peptide_id))
length(as.vector(Lgod_counts_table_ann_up$ensembl_peptide_id))
clade2_common_genes <- union(as.vector(Fparv_counts_table_ann_up$ensembl_peptide_id),as.vector(Lgod_counts_table_ann_up$ensembl_peptide_id))
length(clade2_common_genes)
length(as.vector(Lpar_counts_table_ann_up$ensembl_peptide_id))
clade2_common_genes <- union(clade2_common_genes,as.vector(Lpar_counts_table_ann_up$ensembl_peptide_id))
length(clade2_common_genes)
# Clade 3



# ============================================
#
# Heatmap of common genes
#
# ============================================
# separate between up and down

counts_table = counts(dds, normalized=TRUE)
#clade1_cols <- c(Fcat_cols,Fdia_cols,Fgran_cols,FhetMDPL_cols,FhetMDPP_cols,Frath_cols,Fsim_cols)
Fcat_TR_cols <- c("F_catanatus_transfer_1.quant","F_catanatus_transfer_2.quant")
Fcat_FW_cols <- c("F_catanatus_FW_1.quant","F_catanatus_FW_2.quant")
Fcat_BW_cols <- c("F_catanatus_BW_1.quant","F_catanatus_BW_2.quant","F_catanatus_BW_3.quant")
Fdia_FW_cols <- c("F_diaphanus_FW_2.quant", "F_diaphanus_FW_3.quant")
Fdia_TR_cols <- c("F_diaphanus_transfer_1.quant","F_diaphanus_transfer_2.quant")
Fdia_BW_cols <- c("F_diaphanus_BW_1.quant", "F_diaphanus_BW_2.quant")
Fgran_FW_cols <- c("F_grandis_FW_1.quant", "F_grandis_FW_2.quant", "F_grandis_FW_3.quant")
Fgran_TR_cols <- c("F_grandis_transfer_1.quant", "F_grandis_transfer_2.quant", "F_grandis_transfer_3.quant")
Fgran_BW_cols <- c("F_grandis_BW_1.quant", "F_grandis_BW_2.quant", "F_grandis_BW_3.quant")
FhetMDPL_FW_cols <- c("F_heteroclitusMDPL_FW_1.quant", "F_heteroclitusMDPL_FW_2.quant", "F_heteroclitusMDPL_FW_3.quant")
FhetMDPL_TR_cols <- c("F_heteroclitusMDPL_transfer_1.quant", "F_heteroclitusMDPL_transfer_2.quant", "F_heteroclitusMDPL_transfer_3.quant")
FhetMDPL_BW_cols <- c("F_heteroclitusMDPL_BW_1.quant" , "F_heteroclitusMDPL_BW_2.quant" , "F_heteroclitusMDPL_BW_3.quant")
FhetMDPP_BW_cols <- c("F_heteroclitusMDPP_BW_1.quant" , "F_heteroclitusMDPP_BW_2.quant" , "F_heteroclitusMDPP_BW_3.quant")
FhetMDPP_FW_cols <- c("F_heteroclitusMDPP_FW_1.quant", "F_heteroclitusMDPP_FW_2.quant", "F_heteroclitusMDPP_FW_3.quant")
FhetMDPP_TR_cols <- c("F_heteroclitusMDPP_transfer_1.quant", "F_heteroclitusMDPP_transfer_2.quant", "F_heteroclitusMDPP_transfer_3.quant")
Frath_BW_cols <- c("F_rathbuni_BW_1.quant","F_rathbuni_BW_2.quant","F_rathbuni_BW_3.quant")
Frath_FW_cols <- c("F_rathbuni_FW_1.quant","F_rathbuni_FW_2.quant","F_rathbuni_FW_3.quant")
Frath_TR_cols <- c("F_rathbuni_transfer_1.quant","F_rathbuni_transfer_2.quant","F_rathbuni_transfer_3.quant")
Fsim_BW_cols <- c("F_similis_BW_1.quant","F_similis_BW_2.quant","F_similis_BW_3.quant")
Fsim_FW_cols <- c("F_similis_FW_1.quant","F_similis_FW_2.quant","F_similis_FW_3.quant")
Fsim_TR_cols <- c("F_similis_transfer_1.quant","F_similis_transfer_2.quant","F_similis_transfer_3.quant")

clade1_BW <- c(Fsim_BW_cols,Fgran_BW_cols,FhetMDPL_BW_cols,FhetMDPP_BW_cols,Fdia_BW_cols,Frath_BW_cols,Fcat_BW_cols)
clade1_FW <- c(Fsim_FW_cols,Fgran_FW_cols,FhetMDPL_FW_cols,FhetMDPP_FW_cols,Fdia_FW_cols,Frath_FW_cols,Fcat_FW_cols)

clade1_cols <- c(clade1_FW,clade1_BW)

common_counts <- counts_table[rownames(counts_table) %in% clade1_common_genes,]
clade1_counts <- common_counts[,colnames(common_counts) %in% clade1_cols]
clade1_counts <- clade1_counts[,clade1_cols]
id <-rownames(clade1_counts)
d<-as.matrix(clade1_counts)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="Clade1, FC > 0.5, common sig genes (padj<0.05)",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)

                                  