library(UpSetR)

dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)
dim(sig_threeway)
dim(sig_main_clade)
dim(sig_main_physiology)
dim(sig_main_salinity)


# get species and gene annotation numbers
# (Turns out this wasn't necessary.)

dir <- "/Users/johnsolk/Documents/UCDavis/Whitehead/kfish_expression_July2019/tximport_counts_ensembl/"
files_list <- list.files(dir)
files_list
for (species_file in files_list){
  species <- paste(strsplit(species_file, "_")[[1]][1],strsplit(species_file, "_")[[1]][2],sep="_")
  print(species)
  genes <- read.csv(paste(dir,species_file,sep=""),sep="\t")
  print(length(genes$X))
}

# Dataset
# This is the part that works.

pa <- read.csv("/Users/johnsolk/Documents/UCDavis/Whitehead/kfish_expression_July2019/presence_absence.csv")
pa <- pa[,-c(1)]


upset(pa, sets = c("A_xenica","F_catanatus","F_chrysotus","F_diaphanus","F_grandis","F_heteroclitusMDPL","F_heteroclitusMDPP","F_notatus","F_nottii","F_olivaceous","F_parvapinis","F_rathbuni","F_sciadicus","F_similis","F_zebrinus","L_goodei","L_parva"),mb.ratio = c(0.55, 0.45))

