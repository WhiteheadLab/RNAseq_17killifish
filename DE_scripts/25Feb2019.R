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
filter <- rownames(counts[rowSums(counts >= 0.1) >= 18,])
filter <- rownames(counts[rowSums(counts >= 0.1) >= 100,])
#> dds <- DESeq(dds, full = m1, betaPrior=FALSE)
#using supplied model matrix
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#391 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit #argument with nbinomWaldTest
filter <- rownames(counts[rowSums(counts >= 1) >= 122,])
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

cols<-colnames(filtered_counts)
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

# This works, but has problems converging (data need filtering?):
m1 <- model.matrix(~species + species:condition,ExpDesign)
# This breaks:
#m1 <- model.matrix(~ species + condition + species:condition,ExpDesign)
# this could possibly work (rank = 27, dim(m1) = 34)
#m1 <- model.matrix(~ physiology + clade + condition + physiology:clade:condition,ExpDesign)
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
#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#dds <- nbinomWaldTest(dds, maxit=500)

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
counts_table_filtered9k <- counts(dds, normalized=TRUE)
write.csv(counts_table,"../../Ensembl_counts_normalized_filtered9k_25Feb2019.csv")
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

