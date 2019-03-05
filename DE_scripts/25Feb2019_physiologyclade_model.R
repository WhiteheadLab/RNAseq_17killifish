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
drops <- c("X","Ensembl",
           "F_zebrinus_BW_1.quant","F_zebrinus_BW_2.quant",
           "F_zebrinus_FW_1.quant","F_zebrinus_FW_2.quant",
           "F_notti_FW_1.quant","F_notti_FW_2.quant")
counts<-counts_design[!counts_design$Ensembl == 'Empty',]
rownames(counts)<-counts$Ensembl
design <- design[ , !(names(design) %in% drops)]
counts <- counts[ , !(names(counts) %in% drops)]
dim(design)
dim(counts)

#------------------------

# Filter, all 128 samples have a count of at least 0.01

#------------------------
#filter <- rownames(counts[rowSums(counts >= 0.01) >= 64,])
# thousands of genes not converging
filter <- rownames(counts[rowSums(counts >= 0.01) >= 100,])
#> dds <- DESeq(dds, full = m1, betaPrior=FALSE)
#using supplied model matrix
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#528 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
filter <- rownames(counts[rowSums(counts >= 0.01) >= 122,])
#> dim(counts_table)
#[1] 6500  122
#> dds <- DESeq(dds, full = m1, betaPrior=FALSE)
#using supplied model matrix
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#fitting model and testing
#1 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#but this is missing all but 4 of our genes of interest
filtered_counts <- counts[filter,]
dim(filtered_counts)

#------------------------

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

#------------------------

# Make Design

#------------------------

cols<-colnames(counts)
ExpDesign <- data.frame(row.names=cols,
                        condition=condition, 
                        species = species,
                        physiology = physiology,
                        clade=clade,
                        species_condition=species_condition)
ExpDesign

# suggestions from Amanda
# ~ species + clade + condition
# ~  clade + species_condition
# ~ species + clade + condition + species:condition

m1 <- model.matrix(~ species + condition + clade + species:condition,ExpDesign)

# This works:
#m1 <- model.matrix(~clade + physiology + clade:physiology ,ExpDesign)
# this happens:
#> dds <- DESeq(dds, full = m1, betaPrior=FALSE)
#using supplied model matrix
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#22 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#-- replacing outliers and refitting for 307 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 
#-- original counts are preserved in counts(dds)
#estimating dispersions
#fitting model and testing

colnames(m1)
# run these to remove columns with all 0s
#all.zero <- apply(m1, 2, function(x) all(x==0))
#idx <- which(all.zero)
#m1 <- m1[,-idx]

#use this to check:
Matrix::rankMatrix( m1 )
dim(m1)

all(rownames(ExpDesign) == colnames(filtered_counts))
counts_round<- round(data.matrix(filtered_counts),digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = m1)

# try this, suggested from this:
# https://support.bioconductor.org/p/65091/
#dds <- estimateSizeFactors(dds)
#nc <- counts(dds, normalized=TRUE)
#filtered <- rowSums(nc >= 10) >= 2
#dds <- dds[filter,]

# or try this:
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=1000)

# Error in checkFullRank(modelMatrix) : 
#the model matrix is not full rank, so the model cannot be fit as specified.
#One or more variables or interaction terms in the design formula are linear
#combinations of the others and must be removed.
#Please read the vignette section 'Model matrix not full rank':
#  vignette('DESeq2')

dds <- DESeq(dds, full = m1, betaPrior=FALSE)
# This can take ~50 min to run on Jetstream instance s1.xlarge (CPU: 24, Mem: 60 GB, Disk: 240 GB, Disk: 240 GB root)
# > 1 hr on Macbook Pro, 16 GB RAM

# removes rows that do not converge
ddsClean <- dds[which(mcols(dds)$betaConv),]

dds<-ddsClean
counts_table <- counts(dds, normalized=TRUE)
write.csv(counts_table_filtered122_speciescondition,"../../Ensembl_counts_normalized_speciescondition_filtered122_25Feb2019.csv")

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

pdf("../../multi-ggplot2-catalog_salinity_25Feb2019.pdf",paper="USr",width=13.5, height=8)
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
