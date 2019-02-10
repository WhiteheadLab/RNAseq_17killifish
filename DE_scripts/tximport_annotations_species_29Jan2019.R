library(DESeq2)
library(RColorBrewer)
library(gplots)
library(tximport)

setwd("~/Documents/UCDavis/Whitehead/kfish_salmon")
dir="~/Documents/UCDavis/Whitehead/kfish_salmon/kfish_salmon_by_species/"
files_list = list.files("~/Documents/UCDavis/Whitehead/kfish_salmon/kfish_salmon_by_species/")
files_list
#test one species out first
species = "A_xenica"
species_files = paste(dir,species,sep="")
species_files_list = list.files(species_files)
print(species_files)
print(species_files_list)
files <- file.path(dir,species,species_files_list,"quant.sf")
files
print(file.exists(files))
gene_names_file_ncbi = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene_ncbi/",species,sep="")
gene_names_ncbi <- read.csv(paste(gene_names_file_ncbi,".ncbi.tx2gene.csv",sep=""))
print(dim(gene_names_ncbi))
gene_names_id_ncbi <- gene_names_ncbi[,c(2,3)]
tx2gene_ncbi <- gene_names_id_ncbi
colnames(tx2gene_ncbi)<-c("Name","NCBI")
head(tx2gene_ncbi)
dim(tx2gene_ncbi)
tx2gene_ncbi <- tx2gene_ncbi[!duplicated(tx2gene_ncbi$Name), ]
#tx2gene_ncbi <- tx2gene_ncbi[!duplicated(tx2gene$OG), ]
txi_ncbi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_ncbi)
print(dim(txi_ncbi.salmon$counts))
counts <- txi_ncbi.salmon$counts
colnames(counts) <- species_files_list
head(counts)
counts_matrix = paste(dir,"tximport_counts",species,sep="")
counts_matrix = paste(counts_matrix,"_counts.csv",sep="")
write.csv(counts,counts_matrix)


# all species ncbi
for (species in files_list){
  species_files = paste(dir,species,sep="")
  species_files_list = list.files(species_files)
  print(species_files)
  print(species_files_list)
  files <- file.path(dir,species,species_files_list,"quant.sf")
  files
  print(file.exists(files))
  gene_names_file_ncbi = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene_ncbi/",species,sep="")
  gene_names_ncbi <- read.csv(paste(gene_names_file_ncbi,".ncbi.tx2gene.csv",sep=""))
  print(dim(gene_names_ncbi))
  gene_names_id_ncbi <- gene_names_ncbi[,c(2,3)]
  tx2gene_ncbi <- gene_names_id_ncbi
  colnames(tx2gene_ncbi)<-c("Name","NCBI")
  head(tx2gene_ncbi)
  tx2gene_ncbi <- tx2gene_ncbi[!duplicated(tx2gene_ncbi$Name), ]
  txi_ncbi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_ncbi)
  print(dim(txi_ncbi.salmon$counts))
  counts <- txi_ncbi.salmon$counts
  colnames(counts) <- species_files_list
  counts_matrix = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tximport_counts_ncbi/",species,sep="")
  counts_matrix = paste(counts_matrix,"_counts_ncbi.csv",sep="")
  write.table(counts, file=counts_matrix, quote=FALSE, sep='\t', col.names = NA)
}

# all species evigene
for (species in files_list){
  species_files = paste(dir,species,sep="")
  species_files_list = list.files(species_files)
  print(species_files)
  print(species_files_list)
  files <- file.path(dir,species,species_files_list,"quant.sf")
  files
  print(file.exists(files))
  gene_names_file_ncbi = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene_evigene/",species,sep="")
  gene_names_evigene <- read.csv(paste(gene_names_file_ncbi,".evigene.tx2gene.csv",sep=""))
  print(dim(gene_names_evigene))
  gene_names_id_evigene <- gene_names_evigene[,c(2,3)]
  tx2gene_evigene <- gene_names_id_evigene
  colnames(tx2gene_evigene)<-c("Name","Evigene")
  head(tx2gene_evigene)
  tx2gene_evigene <- tx2gene_evigene[!duplicated(tx2gene_evigene$Name), ]
  txi_evigene.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_evigene)
  print(dim(txi_evigene.salmon$counts))
  counts <- txi_evigene.salmon$counts
  colnames(counts) <- species_files_list
  counts_matrix = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tximport_counts_evigene/",species,sep="")
  counts_matrix = paste(counts_matrix,"_counts_evigene.csv",sep="")
  write.table(counts, file=counts_matrix, quote=FALSE, sep='\t', col.names = NA)
}

# all species ensembl
for (species in files_list){
  species_files = paste(dir,species,sep="")
  species_files_list = list.files(species_files)
  print(species_files)
  print(species_files_list)
  files <- file.path(dir,species,species_files_list,"quant.sf")
  files
  print(file.exists(files))
  gene_names_file_ensembl = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene_ensembl/",species,sep="")
  gene_names_ensembl <- read.csv(paste(gene_names_file_ensembl,".ensembl.tx2gene.csv",sep=""))
  print(dim(gene_names_ensembl))
  gene_names_id_ensembl <- gene_names_ensembl[,c(2,3)]
  tx2gene_ensembl <- gene_names_id_ensembl
  colnames(tx2gene_ensembl)<-c("Name","Ensembl")
  head(tx2gene_ensembl)
  tx2gene_ensembl <- tx2gene_ensembl[!duplicated(tx2gene_ensembl$Name), ]
  txi_ensembl.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene_ensembl)
  print(dim(txi_ensembl.salmon$counts))
  counts <- txi_ensembl.salmon$counts
  colnames(counts) <- species_files_list
  counts_matrix = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tximport_counts_ensembl/",species,sep="")
  counts_matrix = paste(counts_matrix,"_counts_ensembl.csv",sep="")
  write.table(counts, file=counts_matrix, quote=FALSE, sep='\t', col.names = NA)
}