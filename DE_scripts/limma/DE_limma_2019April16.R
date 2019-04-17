library(limma)
library('edgeR')
#setwd("DE_scripts/limma")

# This is the counts with Experimental Design Info in the last 5 rows

if(!file.exists('Ensembl_species_counts_designfactors.csv')){
  download.file("https://osf.io/7vp38/download",'Ensembl_species_counts_designfactors.csv')
}

if(!file.exists('nonzero_clade_physiology_counts_design.csv')){
  download.file("https://osf.io/be7ny/download",'nonzero_clade_physiology_counts_design.csv')
}


if(!file.exists('greater5counts_clade_physiology_counts_design.csv')){
  download.file("https://osf.io/be7ny/download",'greater5counts_clade_physiology_counts_design.csv')
}

counts_design <- read.csv("Ensembl_species_counts_designfactors.csv",stringsAsFactors = TRUE)

# -----------------------
# Format design and counts matrix
# -----------------------

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
#[1] 30466    81

# --------------------
# design cateogories
# --------------------

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
design = model.matrix( ~ physiology + condition + physiology:condition, ExpDesign)

#design1 = model.matrix( ~ -1 + physiology:condition, ExpDesign)

colnames(design)
# check rank of matrix
Matrix::rankMatrix( design )
dim(design)

# ---------------

# normalization

# ---------------

counts_round <- round(data.matrix(counts), digits=0)
counts_round <- head(counts_round, n = 19000)
dim(counts_round)
plot(colSums(t(counts_round)) )

cpm <- cpm(counts_round)
lcpm <- cpm(counts_round, log = TRUE)

table(rowMeans(counts_round) < 300)
plotDensities(lcpm, legend = FALSE, main = "Before filtering")
keep.exprs <- rowSums(cpm > 20) >= 40
x <- counts[keep.exprs,]
dim(x)
x <- round(data.matrix(x), digits=0)
lcpm <- cpm(x, log=TRUE)
plotDensities(lcpm, legend = FALSE, main = "After filtering")

# ---------------

# DE analysis

# ---------------

genes = DGEList(count = counts_round, group = condition_physiology)
genes = calcNormFactors( genes )
vobj = voom( genes, design, plot=TRUE)



# Coefficients not estimable: physiologyM:condition15_ppt 
# Warning message:
#  Partial NA coefficients for 30466 probe(s) 
corfit <- duplicateCorrelation(vobj,design,block=ExpDesign$clade)
# Coefficients not estimable: physiologyM:condition15_ppt 
#There were 42 warnings (use warnings() to see them)
corfit$consensus
# [1] 0.2376267
fitRan <- lmFit(vobj,design,block=ExpDesign$clade,correlation=corfit$consensus)
# Coefficients not estimable: physiologyM:condition15_ppt 
# Warning message:
#  Partial NA coefficients for 30466 probe(s) 
fitRan <- eBayes(fitRan)
#Warning message:
#  In .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  :
#               Estimation of var.prior failed - set to default value
topTable(fitRan,coef=ncol(design))

