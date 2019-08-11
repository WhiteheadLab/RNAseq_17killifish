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

# ============================================

# Get biomart annotations from Ensembl
# https://uswest.ensembl.org/Fundulus_heteroclitus/Info/Index

# ============================================

counts_design <- read.csv('~/Documents/UCDavis/Whitehead/kfish_expression_July2019/Ensembl_species_counts_designfactors.csv',stringsAsFactors = FALSE)

dim(counts_design)
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
query<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','go_id','description'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
head(query)
dim(query)
length(unique(query$ensembl_peptide_id))


# GO only - SALINITY
sal <- query[query$ensembl_peptide_id %in% cons_salinity$ID,]
salGO <- sal[,c('ensembl_peptide_id','go_id')]

sal_up <- query[query$ensembl_peptide_id %in% group1up$ID,]
salGO_up <- sal_up[,c('ensembl_peptide_id','go_id')]

sal_down <- query[query$ensembl_peptide_id %in% group2down$ID,]
salGO_down <- sal_down[,c('ensembl_peptide_id','go_id')]

colnames(salGO_up) <- c("ensembl_peptide_id","go_id")
colnames(salGO_down) <- c("ensembl_peptide_id","go_id")
colnames(salGO) <- c("ensembl_peptide_id","go_id")

dim(salGO_up)
dim(salGO_down)
dim(salGO)

salGO <- salGO[salGO$go_id != "",]
salGO_up <- salGO_up[salGO_up$go_id != "",]
salGO_down <- salGO_down[salGO_down$go_id != "",]

dim(salGO)
dim(salGO_up)
dim(salGO_down)


# GO only - SALINITY-PHYSIOLOGY interacction
salphys <- query[query$ensembl_peptide_id %in% salinity_phys$ID,]
salphysGO <- salphys[,c('ensembl_peptide_id','go_id')]

salphys_up <- query[query$ensembl_peptide_id %in% group1up$ID,]
salphysGO_up <- salphys_up[,c('ensembl_peptide_id','go_id')]

salphys_down <- query[query$ensembl_peptide_id %in% group2down$ID,]
salphysGO_down <- salphys_down[,c('ensembl_peptide_id','go_id')]

colnames(salphysGO_up) <- c("ensembl_peptide_id","go_id")
colnames(salphysGO_down) <- c("ensembl_peptide_id","go_id")
colnames(salphysGO) <- c("ensembl_peptide_id","go_id")

dim(salphysGO_up)
dim(salphysGO_down)
dim(salphysGO)

salphysGO <- salphysGO[salphysGO$go_id != "",]
salphysGO_up <- salphysGO_up[salphysGO_up$go_id != "",]
salphysGO_down <- salphysGO_down[salphysGO_down$go_id != "",]

dim(salphysGO)
dim(salphysGO_up)
dim(salphysGO_down)






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
colnames(TERM2GENE) <- c("GO_ID","gene")
TERM2GENE$gene <- as.character(TERM2GENE$gene)
dim(TERM2GENE)
universe <- universe[universe %in% TERM2GENE$gene]
length(universe)
#df = NULL
#for (GO_ID in TERM2GENE$GO_ID){
#  GOterm <- Term(as.list(GOTERM[GO_ID])[[1]])
#  df = rbind(df, data.frame(GO_ID,GOterm))
#}
dim(df)
head(df)

# change gene for each effect
salGO_up$ensembl_peptide_id <- as.character(salGO_up$ensembl_peptide_id)
gene <- salGO_up$ensembl_peptide_id[salGO_up$ensembl_peptide_id %in% universe]

salGO_down$ensembl_peptide_id <- as.character(salGO_down$ensembl_peptide_id)
gene <- salGO_down$ensembl_peptide_id[salGO_down$ensembl_peptide_id %in% universe]

salGO$ensembl_peptide_id <- as.character(salGO$ensembl_peptide_id)
gene <- salGO$ensembl_peptide_id[salGO$ensembl_peptide_id %in% universe]

# PHYSIOLOGY-SALINITY response
salphysGO_up$ensembl_peptide_id <- as.character(salphysGO_up$ensembl_peptide_id)
gene <- salphysGO_up$ensembl_peptide_id[salphysGO_up$ensembl_peptide_id %in% universe]

salphysGO_down$ensembl_peptide_id <- as.character(salphysGO_down$ensembl_peptide_id)
gene <- salphysGO_down$ensembl_peptide_id[salphysGO_down$ensembl_peptide_id %in% universe]

salphysGO$ensembl_peptide_id <- as.character(salphysGO$ensembl_peptide_id)
gene <- salphysGO$ensembl_peptide_id[salphysGO$ensembl_peptide_id %in% universe]


gene<-as.character(unique(gene))
length(gene)
output<-enricher(gene, pvalueCutoff = 0.2, pAdjustMethod = "BH", universe,
                 minGSSize = 2, maxGSSize = 500, qvalueCutoff = 0.2, 
                 TERM2GENE=TERM2GENE,
                 TERM2NAME = df)
dotplot(output)












# find all circadian regulation of gene expression genes
# df has the GOterm
# which GO_ID corresponds to "circadian regulation of gene expression"
circ_go <- df[df$GOterm=="circadian regulation of gene expression",]
circ_reg_go <- df[df$GOterm=="regulation of circadian rhythm",]
circ <- as.character(unique(circ_go$GO_ID))
circ_reg <- as.character(unique(circ_reg_go$GO_ID))
# GO:0032922
# which ID have that GO_ID? in salphysGO_down
circ_genes <- salphys_down[salphys_down$go_id == circ,]$ensembl_gene_id
circ_reg_genes <- salphys_down[salphys_down$go_id == circ_reg,]$ensembl_gene_id
circ_genes <- c(circ_genes, circ_reg_genes) 


# edit GOterm in df to make entry shorter


ox <- df[df$GOterm=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",]
ox_gens <- salphys_down[salphys_down$go_id == circ_reg,]$ensembl_gene_id

test <- df
test[test$GOterm == "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",]$GOterm <- "oxidoreductase activity"

test[test$GO_ID =="GO:0016712",]
df<-test
