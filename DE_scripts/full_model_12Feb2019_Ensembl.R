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
packages(dplyr)
packages(tidyr)
packages(ggplot2)
bioconductors(DESeq2)
bioconductors(tximport)
bioconductors(biomaRt)

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


filter <- rownames(counts[rowSums(counts >= 0.01) >= 128,])
filtered_counts <- counts[filter,]
dim(filtered_counts)

#------------------------

# Make Design categories

#------------------------

# design cateogories (full)
species<-as.character(unlist(design[1,]))
nativephysiology<-as.character(unlist(design[2,]))
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
ExpDesign <- data.frame(row.names=cols, condition=condition, species = species, clade=clade,species_condition=species_condition,species_clade = species_clade)
ExpDesign


m1 <- model.matrix(~condition + species_clade + species_clade:condition,ExpDesign)
colnames(m1)
all.zero <- apply(m1, 2, function(x) all(x==0))
idx <- which(all.zero)
m1 <- m1[,-idx]

all(rownames(ExpDesign) == colnames(filtered_counts))
counts_round<- round(data.matrix(filtered_counts),digits=0)
dds <- DESeqDataSetFromMatrix(countData = counts_round,colData = ExpDesign,design = m1)
dds <- DESeq(dds, full = m1, betaPrior=FALSE)
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
plotDispEsts(dds)

plotPCA(vsd, intgroup=c("species_clade"))
plotPCAWithSampleNames(vsd,intgroup=c("species_clade"),ntop=40000)

# ---------------------------
# make a matrix of counts
# ---------------------------

counts_table = counts(dds, normalized=TRUE)

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
       
      
    
        
  



                                  