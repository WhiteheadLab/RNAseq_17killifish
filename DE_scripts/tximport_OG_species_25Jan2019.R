library(DESeq2)
library(RColorBrewer)
library(gplots)
library(tximport)

setwd("~/Documents/UCDavis/Whitehead/kfish_salmon")
dir="~/Documents/UCDavis/Whitehead/kfish_salmon/salmon_renamed_species/"
files_list = list.files("~/Documents/UCDavis/Whitehead/kfish_salmon/salmon_renamed_species/")
files_list
#test one species out first
species = "L_parva"
species_files = paste(dir,species,sep="")
species_files_list = list.files(species_files)
print(species_files)
print(species_files_list)
files <- file.path(dir,species,species_files_list)
files
print(file.exists(files))
gene_names_file = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene/",species,sep="")
gene_names <- read.csv(paste(gene_names_file,"_tx2gene.csv",sep=""))
print(dim(gene_names))
gene_names_id <- gene_names[,c(2,3)]
tx2gene <- gene_names_id
colnames(tx2gene)<-c("Name","OG")
head(tx2gene)
tx2gene <- tx2gene[!duplicated(tx2gene$transcriptID), ]
tx2gene <- tx2gene[!duplicated(tx2gene$OG), ]
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
print(dim(txi.salmon$counts))
counts <- txi.salmon$counts
colnames(counts) <- species_files_list
head(counts)
counts_matrix = paste(dir,"tximport_counts",species,sep="")
counts_matrix = paste(counts_matrix,"_counts.csv",sep="")
write.csv(counts,counts_matrix)

test_csv<-readr::read_tsv("~/Documents/UCDavis/Whitehead/kfish_salmon/salmon_renamed_species//A_xenica/A_xenica_transfer_3.quant.sf",progress=FALSE, col_types=col.types)
head(test_csv)



# all species
for (species in files_list){
  species_files = paste(dir,species,sep="")
  species_files_list = list.files(species_files)
  print(species_files)
  print(species_files_list)
  files <- file.path(dir,species,species_files_list)
  files
  print(file.exists(files))
  gene_names_file = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tx2gene/",species,sep="")
  gene_names <- read.csv(paste(gene_names_file,"_tx2gene.csv",sep=""))
  print(dim(gene_names))
  gene_names_id <- gene_names[,c(2,3)]
  tx2gene <- gene_names_id
  tx2gene <- tx2gene[!duplicated(tx2gene$transcriptID), ]
  tx2gene <- tx2gene[!duplicated(tx2gene$OG), ]
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
  print(dim(txi.salmon$counts))
  counts <- txi.salmon$counts
  colnames(counts) <- species_files_list
  head(counts)
  counts_matrix = paste("~/Documents/UCDavis/Whitehead/kfish_salmon/tximport_counts/",species,sep="")
  counts_matrix = paste(counts_matrix,"_counts.csv",sep="")
  #write.csv(counts,counts_matrix)
  write.table(counts, file=counts_matrix, quote=FALSE, sep='\t', col.names = NA)
}
