#------------------------

# Load packages

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
library(goseq)
library(GO.db)
library(AnnotationDbi)
library(clusterProfiler)
setwd("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_scripts")
# ------------------------

# Import data from each species DE analysis

# ------------------------
# load each species response table
# ------------------------
A_xenica_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/A_xenica_stats_annotations_counts.csv")
F_catenatus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_catanatus_stats_annotations_counts.csv")
F_chrysotus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_chrysotus_stats_annotations_counts.csv")
F_diaphanus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_diaphanus_stats_annotations_counts.csv")
F_grandis_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_grandis_stats_annotations_counts.csv")
F_heteroclitusMDPL_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_heteroclitusMDPL_stats_annotations_counts.csv")
F_heteroclitusMDPP_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_heteroclitusMDPP_stats_annotations_counts.csv")
F_notatus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_notatus_stats_annotations_counts.csv")
F_olivaceus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_olivaceous_stats_annotations_counts.csv")
F_parvapinis_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_parvapinis_stats_annotations_counts.csv")
F_rathbuni_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_rathbuni_stats_annotations_counts.csv")
F_sciadicus_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_sciadicus_stats_annotations_counts.csv")
F_similis_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/F_similis_stats_annotations_counts.csv")
L_goodei_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/L_goodei_stats_annotations_counts.csv")
L_parva_stats <- read.csv("~/Documents/UCDavis/Whitehead/counts_stats_byspecies/L_parva_stats_annotations_counts.csv")
# --------------------------------
# isolate only padj<0.05 for each
# --------------------------------
A_xenica_stats_sig <- subset(A_xenica_stats, A_xenica_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(A_xenica_stats_sig)
F_catenatus_stats_sig <- subset(F_catenatus_stats, F_catenatus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_catenatus_stats_sig)
F_chrysotus_stats_sig <- subset(F_chrysotus_stats, F_chrysotus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_chrysotus_stats_sig)
F_diaphanus_stats_sig <- subset(F_diaphanus_stats, F_diaphanus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_diaphanus_stats_sig)
F_grandis_stats_sig <- subset(F_grandis_stats, F_grandis_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_grandis_stats_sig)
F_heteroclitusMDPL_stats_sig <- subset(F_heteroclitusMDPL_stats, F_heteroclitusMDPL_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_heteroclitusMDPL_stats_sig)
F_heteroclitusMDPP_stats_sig <- subset(F_heteroclitusMDPP_stats, F_heteroclitusMDPP_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_heteroclitusMDPP_stats_sig)
F_notatus_stats_sig <- subset(F_notatus_stats, F_notatus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_notatus_stats_sig)
F_olivaceus_stats_sig <- subset(F_olivaceus_stats, F_olivaceus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_olivaceus_stats_sig)
F_parvapinis_stats_sig <- subset(F_parvapinis_stats, F_parvapinis_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_parvapinis_stats_sig)
F_rathbuni_stats_sig <- subset(F_rathbuni_stats, F_rathbuni_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_rathbuni_stats_sig)
F_sciadicus_stats_sig <- subset(F_sciadicus_stats, F_sciadicus_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_sciadicus_stats_sig)
F_similis_stats_sig <- subset(F_similis_stats, F_similis_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(F_similis_stats_sig)
L_goodei_stats_sig <- subset(L_goodei_stats, L_goodei_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(L_goodei_stats_sig)
L_parva_stats_sig <- subset(L_parva_stats, L_parva_stats$`padj.15ppt.v.0.2ppt`<= 0.05)
dim(L_parva_stats_sig)

# ============================================

# Get biomart annotations from Ensembl
# https://uswest.ensembl.org/Fundulus_heteroclitus/Info/Index

# ============================================

counts_design <- read.csv("../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)
design <- counts_design[counts_design$Ensembl == 'Empty',]
#design$type <- c("species","native_salinity","clade","group","condition")
drops <- c("X","Ensembl")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
dim(design)
dim(counts)

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
ensembl_proteinID = rownames(counts)
length(ensembl_proteinID)
query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','go_id','description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
head(query)
dim(query)
length(unique(query$ensembl_peptide_id))

# ============================================
#
# Get GO terms
# combine each sig table with annotations
# ============================================

AxenGO_all <- merge(A_xenica_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FcatGO_all <- merge(F_catenatus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FchryGO_all <- merge(F_chrysotus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FdiaGO_all <- merge(F_diaphanus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FgranGO_all <- merge(F_grandis_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FhetMDPLGO_all <- merge(F_heteroclitusMDPL_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FhetMDPPGO_all <- merge(F_heteroclitusMDPP_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FnotaGO_all <- merge(F_notatus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FoliGO_all <- merge(F_olivaceus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FparvGO_all <- merge(F_parvapinis_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FrathGO_all <- merge(F_rathbuni_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FsciGO_all <- merge(F_sciadicus_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
FsimGO_all <- merge(F_similis_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
LgodGO_all <- merge(L_goodei_stats_sig,query,by.x='X',by.y='ensembl_peptide_id') 
LparGO_all <- merge(L_parva_stats_sig,query,by.x='X',by.y='ensembl_peptide_id')
# GO only
AxenGO <- AxenGO_all[,c('go_id', 'ensembl_peptide_id')]
FcatGO <- FcatGO_all[,c('go_id', 'ensembl_peptide_id')]
FchryGO <- FchryGO_all[,c('go_id', 'ensembl_peptide_id')]
FdiaGO <- FdiaGO_all[,c('go_id', 'ensembl_peptide_id')]
FgranGO <- FgranGO_all[,c('go_id', 'ensembl_peptide_id')]
FhetMDPLGO <- FhetMDPLGO_all[,c('go_id', 'ensembl_peptide_id')]
FhetMDPPGO <- FhetMDPPGO_all[,c('go_id', 'ensembl_peptide_id')]
FnotaGO <- FnotaGO_all[,c('go_id', 'ensembl_peptide_id')]
FoliGO <- FoliGO_all[,c('go_id', 'ensembl_peptide_id')]
FparvGO <- FparvGO_all[,c('go_id', 'ensembl_peptide_id')]
FrathGO <- FrathGO_all[,c('go_id', 'ensembl_peptide_id')]
FsciGO <- FsciGO_all[,c('go_id', 'ensembl_peptide_id')]
FsimGO <- FsimGO_all[,c('go_id', 'ensembl_peptide_id')]
LgodGO <- LgodGO_all[,c('go_id', 'ensembl_peptide_id')]
LparGO <- LparGO_all[,c('go_id', 'ensembl_peptide_id')]
# gene = character vector with ENSEMBL ID of differentially expressed genes 
# universe = character vector with ENSEMBL ID of any gene in transcriptome
# TERM2GENE = mapping bewteen Ensembl ID and GOID from universe
# TERM2NAME = df with GO ID and name of term from the universe


# THIS WORKS to look up terms with GO ID
# Term(as.list(GOTERM["GO:0000166"])[[1]])

# these don't change
universe_genes<-query[,c('ensembl_peptide_id','go_id')]
dim(universe_genes)
universe_genes <- universe_genes[universe_genes$go_id != "",]
dim(universe_genes)
universe <- as.character(unique(universe_genes$ensembl_peptide_id))
length(universe)

# these don't change
TERM2GENE <- query[,c('go_id','ensembl_peptide_id')]
dim(TERM2GENE)
TERM2GENE <- TERM2GENE[TERM2GENE$go_id != "",]
dim(TERM2GENE)
TERM2GENE$gene <- as.character(TERM2GENE$ensembl_peptide_id)
colnames(TERM2GENE) <- c("GO_ID","gene")
dim(TERM2GENE)
universe <- universe[universe %in% TERM2GENE$gene]
length(universe)
df = NULL
for (GO_ID in TERM2GENE$GO_ID){
  GOterm <- Term(as.list(GOTERM[GO_ID])[[1]])
  df = rbind(df, data.frame(GO_ID,GOterm))
}
dim(df)
head(df)

# change gene for each species
gene <- FcatGO$ensembl_peptide_id[FcatGO$ensembl_peptide_id %in% universe]
gene<-as.character(unique(gene))
length(gene)
output<-enricher(gene, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)
dotplot(output)
#View(output@result)
species_list<-list(AxenGO,FcatGO,FchryGO,FdiaGO,FgranGO,FhetMDPLGO,
                   FhetMDPPGO,FnotaGO,FoliGO,FparvGO,FrathGO,FsciGO,
                   FsimGO,LgodGO,LparGO)
species_names <- c("Axen","Fcat","Fchry","Fdia","Fgran","FhetMDPL",
                   "FhetMDPP","Fnota","Foli","Fparv","Frath","Fsci",
                   "Fsim","Lgod","Lpar")


dot_plots <- list()
for (i in 1:length(species_list)){
  print(head(species_list[[i]]))
  gene <- species_list[[i]]$ensembl_peptide_id[species_list[[i]]$ensembl_peptide_id %in% universe]
  gene<-as.character(unique(gene))
  print(length(gene))
  print(head(gene))
  output<-enricher(gene, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                   minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                   TERM2GENE=TERM2GENE,
                   TERM2NAME = df)
  dot_plots[[i]] <- dotplot(output,title=species_names[i])
}
pdf("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_analysis/dotplots.pdf",paper="USr",width=13.5, height=8)
dot_plots
dev.off()

# separate by clade and physiology
# then combine
# Clade 1, FW
# Fcat
# Frath
cladelist <- list(FcatGO,FrathGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
  }
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)

dotplot(output,title="Clade1, FW")


FWclade1 <- 
# Clade 1, M
# Fsim
# FhetMDPP
# FhetMDPL
# Fgran
# Fdia
cladelist <- list(FsimGO,FhetMDPLGO,FhetMDPPGO,FgranGO,FdiaGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
}
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)
dotplot(output,title="Clade1, M")
View(output@result)

# Clade 2, FW
# Lgod
cladelist <- list(LgodGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
}
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)
dotplot(output,title="Clade2, FW (L. goodei)")
View(output@result)
# Clade 2, M
# Fparv
# Lpar
cladelist <- list(FparvGO,LparGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
}
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)

dotplot(output,title="Clade2, M")

# Clade 3, FW
# Foli
# Fnota
# Fsci
cladelist <- list(FoliGO,FnotaGO,FsciGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
}
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)

dotplot(output,title="Clade3, FW")
View(output@result)
write.csv(output@result,"~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_analysis/GOenrichment_Clade3FW.csv")
# Clade 3, M
# Fchry
# Axen
cladelist <- list(FchryGO,AxenGO)
combined <- c()
for (i in 1:length(cladelist)){
  gene <- cladelist[[i]]$ensembl_peptide_id[cladelist[[i]]$ensembl_peptide_id %in% universe]
  gene <- as.character(unique(gene))
  combined <- c(combined,gene)
  combined <- as.character(unique(combined))
}
length(combined)
output<-enricher(combined, pvalueCutoff = 0.1, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)

dotplot(output,title="Clade3, M")
View(output@result)



