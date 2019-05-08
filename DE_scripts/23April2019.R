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
library(BiocParallel)
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
counts_design <- read.csv("../../../Ensembl_species_counts_designfactors.csv",stringsAsFactors = FALSE)
dim(counts_design)
#[1] 30471   130
#------------------------

# Format design and counts

#------------------------

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

#------------------------

# Make Design categories

#------------------------

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
design = model.matrix( ~ physiology*condition*clade, ExpDesign)

colnames(design)
# check rank of matrix
Matrix::rankMatrix( design )
dim(design)


# ---------------------------
# Run DESeq2
# ---------------------------
counts<-as.matrix(as.data.frame(sapply(counts, as.numeric)))
rownames(counts)<-gene.names
all(rownames(ExpDesign) == colnames(counts))
counts_round<- round(counts,digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = design)

dds <- DESeq(dds, full = design, parallel=TRUE, betaPrior=FALSE)
# took 10 min
#using supplied model matrix
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates: 6 workers
#mean-dispersion relationship
#final dispersion estimates, fitting model and testing: 6 workers
#-- replacing outliers and refitting for 407 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing
#13 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit #argument with nbinomWaldTest
# removes rows that do not converge
ddsClean <- dds[which(mcols(dds)$betaConv),]

dds<-ddsClean
counts_table <- counts(dds, normalized=TRUE)
write.csv(counts_table,"../../Ensembl_DESeq2counts_normalized_23April2019.csv")

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
plotPCA(vsd, intgroup=c("species"))
plotPCAWithSampleNames(vsd,intgroup=c("species"))

# ============================================
#
# Genes of Interest
# This was very helpful:
# https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf
# 
# ============================================

all_goi<-c("ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
           "ENSFHEP00000015383","ENSFHEP00000009753","ENSFHEP00000006725","ENSFHEP00000008393",
           "ENSFHEP00000013324","ENSFHEP00000001609","ENSFHEP00000013324","ENSFHEP00000034177",
           "ENSFHEP00000015765","ENSFHEP00000017303","ENSFHEP00000000036","ENSFHEP00000031108",
           "ENSFHEP00000016853","ENSFHEP00000003908")

pdf("../../multi-ggplot2-catalog_salinity_23April2019.pdf",paper="USr",width=13.5, height=8)
for (i in all_goi){
  if (i %in% rownames(counts_table)) {
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
  else {
    print("Not present in filtered data set:")
    print(i)
  }
}
dev.off()


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

# =============================
# Results, contrasts
# =============================

resultsNames(dds)
res_physiologyM.condition15_ppt.cladeClade3_v_Cl1_15ppt_M <- results(dds, tidy=TRUE, contrast("Intercept","physiologyM.condition15_ppt.cladeClade3")) %>% arrange(padj) %>% tbl_df()

# all others
res_Cl3__v_Cl2_15ppt_M <- results(dds, tidy=TRUE, name="physiologyM.condition15_ppt.cladeClade3") %>% arrange(padj) %>% tbl_df()

