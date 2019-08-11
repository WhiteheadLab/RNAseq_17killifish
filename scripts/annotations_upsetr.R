library(UpSetR)

dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)
dim(sig_threeway)
dim(sig_main_clade)
dim(sig_main_physiology)
dim(sig_main_salinity)

input <- c()

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
input <- c(A_xenica = 21733,
           F_catanatus = 22226,
           F_chrysotus = 22222,
           F_diaphanus = 21457,
           F_grandis = 24245,
           F_heteroclitusMDPL = 22806,
           F_heteroclitusMDPP = 22922,
           F_notatus = 22290,
           F_nottii = 18113,
           F_olivaceous = 21482,
           F_parvapinis = 20432,
           F_rathbuni = 22781,
           F_sciadicus = 19881,
           F_similis = 21828,
           F_zebrinus = 20275,
           L_goodei = 21897,
           L_parva = 22157)

pa <- read.csv("/Users/johnsolk/Documents/UCDavis/Whitehead/kfish_expression_July2019/presence_absence.csv")
pa <- pa[,-c(1)]


upset(pa, sets = c("A_xenica","F_catanatus","F_chrysotus","F_diaphanus","F_grandis","F_heteroclitusMDPL","F_heteroclitusMDPP","F_notatus","F_nottii","F_olivaceous","F_parvapinis","F_rathbuni","F_sciadicus","F_similis","F_zebrinus","L_goodei","L_parva"),mb.ratio = c(0.55, 0.45))

