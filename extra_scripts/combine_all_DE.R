library(RColorBrewer)
library(gplots)
# ALL genes
setwd("~/Documents/UCDavis/Whitehead/DE_results_all_genes")
A_xenica_all = read.csv("A_xenica_DE_all.csv")
A_xenica_data <- A_xenica_all[,c(2:11)]
#A_xenica_DE <- A_xenica_all[A_xenica_all$GeneID %in% DE_genes_FW,]
# 8 samples
F_chrysotus_all = read.csv("F_chrysotus_DE_all.csv")
F_chrysotus_data <- F_chrysotus_all[,c(2:10)]
F_grandis_all = read.csv("F_grandis_DE_all.csv")
F_grandis_data <- F_grandis_all[,c(2:11)]
F_heteroclitusMDPL_all = read.csv("F_heteroclitusMDPL_DE_all.csv")
F_heteroclitusMDPL_data <- F_heteroclitusMDPL_all[,c(2:11)]
F_heteroclitusMDPP_all = read.csv("F_heteroclitusMDPP_DE_all.csv")
F_heteroclitusMDPP_data <- F_heteroclitusMDPP_all[,c(2:11)]
F_notatus_all = read.csv("F_notatus_DE_all.csv")
F_notatus_data <- F_notatus_all[,c(2:11)]
# 8 samples
F_olivaceus_all = read.csv("F_olivaceus_DE_all.csv")
F_olivaceus_data <- F_olivaceus_all[,c(2:10)]
# 8 samples
F_parvapinis_all = read.csv("F_parvapinis_DE_all.csv")
F_parvapinis_data <- F_parvapinis_all[,c(2:10)]
F_rathbuni_all = read.csv("F_rathbuni_DE_all.csv")
F_rathbuni_data <- F_rathbuni_all[,c(2:11)]
F_similis_all = read.csv("F_similis_DE_all.csv")
F_similis_data <- F_similis_all[,c(2:11)]
L_goodei_all = read.csv("L_goodei_DE_all.csv")
L_goodei_data <- L_goodei_all[,c(2:11)]
L_parva_all = read.csv("L_parva_DE_all.csv")
L_parva_data <- L_parva_all[,c(2:11)]
remove(L_goodei_all)
remove(L_parva_all)
remove(F_similis_all)
remove(F_rathbuni_all)
remove(F_olivaceus_all)
remove(F_notatus_all)
remove(F_heteroclitusMDPP_all)
remove(F_heteroclitusMDPL_all)
remove(F_grandis_all)
remove(F_chrysotus_all)
remove(A_xenica_all)
remove(F_parvapinis_all)
# Get GeneID in common
F_similis_GeneID<-F_similis_data$GeneID
A_xenica_GeneID <- A_xenica_data$GeneID
length(unique(A_xenica_GeneID))
F_chrysotus_GeneID <- F_chrysotus_data$GeneID
length(unique(F_chrysotus_GeneID))
F_grandis_GeneID <- F_grandis_data$GeneID
length(unique(F_grandis_GeneID))
F_heteroclitusMDPL_GeneID <- F_heteroclitusMDPL_data$GeneID
length(unique(F_grandis_GeneID))
F_heteroclitusMDPP_GeneID <- F_heteroclitusMDPP_data$GeneID
length(unique(F_heteroclitusMDPP_GeneID))
F_notatus_GeneID <- F_notatus_data$GeneID
length(unique(F_notatus_GeneID))
F_olivaceus_GeneID <- F_olivaceus_data$GeneID
length(unique(F_olivaceus_GeneID))
F_parvapinis_GeneID <- F_parvapinis_data$GeneID
length(unique(F_parvapinis_GeneID))
F_rathbuni_GeneID <- F_rathbuni_data$GeneID
length(unique(F_rathbuni_GeneID))
L_goodei_GeneID <- L_goodei_data$GeneID
length(unique(L_goodei_GeneID))
L_parva_GeneID <- L_parva_data$GeneID
length(unique(L_parva_GeneID))
all_GeneID<-Reduce(intersect, list(F_similis_GeneID,L_parva_GeneID,L_goodei_GeneID,F_rathbuni_GeneID,F_parvapinis_GeneID,F_olivaceus_GeneID,
                                   F_notatus_GeneID,F_heteroclitusMDPP_GeneID,F_heteroclitusMDPL_GeneID,F_grandis_GeneID,
                                   F_chrysotus_GeneID,A_xenica_GeneID))
sub_F_similis_data<-F_similis_data[F_similis_data$GeneID %in% all_GeneID,]
sub_F_similis_data<-sub_F_similis_data[!duplicated(sub_F_similis_data$GeneID),]
sub_L_parva_data<-L_parva_data[L_parva_data$GeneID %in% all_GeneID,]
sub_L_parva_data<-sub_L_parva_data[!duplicated(sub_L_parva_data$GeneID),]
sub_L_goodei_data<-L_goodei_data[L_goodei_data$GeneID %in% all_GeneID,]
sub_L_goodei_data<-sub_L_goodei_data[!duplicated(sub_L_goodei_data$GeneID),]
sub_F_rathbuni_data<-F_rathbuni_data[F_rathbuni_data$GeneID %in% all_GeneID,]
sub_F_rathbuni_data<-sub_F_rathbuni_data[!duplicated(sub_F_rathbuni_data$GeneID),]
sub_F_parvapinis_data<-F_parvapinis_data[F_parvapinis_data$GeneID %in% all_GeneID,]
sub_F_parvapinis_data<-sub_F_parvapinis_data[!duplicated(sub_F_parvapinis_data$GeneID),]
sub_F_olivaceus_data<-F_olivaceus_data[F_olivaceus_data$GeneID %in% all_GeneID,]
sub_F_olivaceus_data<-sub_F_olivaceus_data[!duplicated(sub_F_olivaceus_data$GeneID),]
sub_F_notatus_data<-F_notatus_data[F_notatus_data$GeneID %in% all_GeneID,]
sub_F_notatus_data<-sub_F_notatus_data[!duplicated(sub_F_notatus_data$GeneID),]
sub_F_heteroclitusMDPP_data<-F_heteroclitusMDPP_data[F_heteroclitusMDPP_data$GeneID %in% all_GeneID,]
sub_F_heteroclitusMDPP_data<-sub_F_heteroclitusMDPP_data[!duplicated(sub_F_heteroclitusMDPP_data$GeneID),]
sub_F_heteroclitusMDPL_data<-F_heteroclitusMDPL_data[F_heteroclitusMDPL_data$GeneID %in% all_GeneID,]
sub_F_heteroclitusMDPL_data<-sub_F_heteroclitusMDPL_data[!duplicated(sub_F_heteroclitusMDPL_data$GeneID),]
sub_F_grandis_data<-F_grandis_data[F_grandis_data$GeneID %in% all_GeneID,]
sub_F_grandis_data<-sub_F_grandis_data[!duplicated(sub_F_grandis_data$GeneID),]
sub_F_chrysotus_data<-F_chrysotus_data[F_chrysotus_data$GeneID %in% all_GeneID,]
sub_F_chrysotus_data<-sub_F_chrysotus_data[!duplicated(sub_F_chrysotus_data$GeneID),]
sub_A_xenica_data<-A_xenica_data[A_xenica_data$GeneID %in% all_GeneID,]
sub_A_xenica_data<-sub_A_xenica_data[!duplicated(sub_A_xenica_data$GeneID),]
# Combine all
all <- merge(sub_F_rathbuni_data,sub_F_parvapinis_data,by="GeneID")
all <- merge(all,sub_A_xenica_data,by="GeneID")             
all <- merge(all,sub_F_chrysotus_data,by="GeneID")             
all <- merge(all,sub_F_grandis_data,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPL_data,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPP_data,by="GeneID")             
all <- merge(all,sub_F_notatus_data,by="GeneID")             
all <- merge(all,sub_F_olivaceus_data,by="GeneID")             
all <- merge(all,sub_L_goodei_data,by="GeneID")             
all <- merge(all,sub_L_parva_data,by="GeneID")
all <- merge(all,sub_F_similis_data,by="GeneID")
dim(all)
colnames(all)
# remove GeneID column
all<-all[,c(2:106)]
colnames(all)
# order samples by physiology
# FW: F_notatus,F_olivaceus,F_rathbuni,L_goodei
# BW: F_parvapinis,F_chrysotus,
# M: F_grandis,F_heteroclitusMDPL,F_heteroclitusMDPP,F_similis,A_xenica, L_parva
all_ordered<-all[,c(71:78,62:70,1:9,
                    18:26,27:34,88:96,
                    35:43,44:61,97:105,79:87,10:17)]
d<-as.matrix(all_ordered)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="Fundulus (ordered by physiology), common genes", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)

setwd("~/Documents/UCDavis/Whitehead/RNAseq_15killifish/DE_results")
# Get union of significant DE genes for all species (inclusive of all, regardless of unique or in common)
# FW
A_xenica_FW = read.csv("A_xenica_DE_TR_FW.csv")
A_xenica_DE = unique(A_xenica_FW$GeneID)
length(A_xenica_DE)
F_chrysotus_FW = read.csv("F_chrysotus_DE_TR_FW.csv")
F_chrysotus_DE = unique(F_chrysotus_FW$GeneID)
length(F_chrysotus_DE)
DE_genes = union(A_xenica_DE,F_chrysotus_DE)
length(DE_genes)
F_grandis_FW = read.csv("F_grandis_DE_TR_FW.csv")
F_grandis_DE = unique(F_grandis_FW$GeneID)
length(F_grandis_DE)
DE_genes = union(DE_genes,F_grandis_DE)
length(DE_genes)
F_heteroclitusMDPL_FW = read.csv("F_heteroclitusMDPL_DE_TR_FW.csv")
F_heteroclitusMDPL_DE = unique(F_heteroclitusMDPL_FW$GeneID)
length(F_heteroclitusMDPL_DE)
DE_genes = union(DE_genes,F_heteroclitusMDPL_DE)
length(DE_genes)
F_heteroclitusMDPP_FW = read.csv("F_heteroclitusMDPP_DE_TR_FW.csv")
F_heteroclitusMDPP_DE = unique(F_heteroclitusMDPP_FW$GeneID)
length(F_heteroclitusMDPP_DE)
DE_genes = union(DE_genes,F_heteroclitusMDPP_DE)
length(DE_genes)
F_notatus_FW = read.csv("F_notatus_DE_TR_FW.csv")
F_notatus_DE = unique(F_notatus_FW$GeneID)
length(F_notatus_DE)
DE_genes = union(DE_genes,F_notatus_DE)
length(DE_genes)
F_olivaceus_FW = read.csv("F_olivaceus_DE_TR_FW.csv")
F_olivaceus_DE = unique(F_olivaceus_FW$GeneID)
length(F_olivaceus_DE)
DE_genes = union(DE_genes,F_olivaceus_DE)
length(DE_genes)
F_parvapinis_FW = read.csv("F_parvapinis_DE_TR_FW.csv")
F_parvapinis_DE = unique(F_parvapinis_FW$GeneID)
length(F_parvapinis_DE)
DE_genes = union(DE_genes,F_parvapinis_DE)
length(DE_genes)
F_rathbuni_FW = read.csv("F_rathbuni_DE_TR_FW.csv")
F_rathbuni_DE = unique(F_rathbuni_FW$GeneID)
length(F_rathbuni_DE)
DE_genes = union(DE_genes,F_rathbuni_DE)
length(DE_genes)
F_similis_FW = read.csv("F_similis_DE_TR_FW.csv")
F_similis_DE = unique(F_similis_FW$GeneID)
length(F_similis_DE)
DE_genes = union(DE_genes,F_similis_DE)
length(DE_genes)
L_goodei_FW = read.csv("L_goodei_DE_TR_FW.csv")
L_goodei_DE = unique(L_goodei_FW$GeneID)
length(L_goodei_DE)
DE_genes = union(DE_genes,L_goodei_DE)
length(DE_genes)
L_parva_FW = read.csv("L_parva_DE_TR_FW.csv")
L_parva_DE = unique(L_parva_FW$GeneID)
length(L_parva_DE)
DE_genes = union(DE_genes,L_parva_DE)
length(DE_genes)
DE_genes_FW <- DE_genes
# BW
A_xenica_BW = read.csv("A_xenica_DE_TR_BW.csv")
A_xenica_DE = unique(A_xenica_BW$GeneID)
length(A_xenica_DE)
F_chrysotus_BW = read.csv("F_chrysotus_DE_TR_BW.csv")
F_chrysotus_DE = unique(F_chrysotus_BW$GeneID)
length(F_chrysotus_DE)
DE_genes = union(A_xenica_DE,F_chrysotus_DE)
length(DE_genes)
F_grandis_BW = read.csv("F_grandis_DE_TR_FW.csv")
F_grandis_DE = unique(F_grandis_BW$GeneID)
length(F_grandis_DE)
DE_genes = union(DE_genes,F_grandis_DE)
length(DE_genes)
F_heteroclitusMDPL_BW = read.csv("F_heteroclitusMDPL_DE_TR_BW.csv")
F_heteroclitusMDPL_DE = unique(F_heteroclitusMDPL_BW$GeneID)
length(F_heteroclitusMDPL_DE)
DE_genes = union(DE_genes,F_heteroclitusMDPL_DE)
length(DE_genes)
F_heteroclitusMDPP_BW = read.csv("F_heteroclitusMDPP_DE_TR_BW.csv")
F_heteroclitusMDPP_DE = unique(F_heteroclitusMDPP_BW$GeneID)
length(F_heteroclitusMDPP_DE)
DE_genes = union(DE_genes,F_heteroclitusMDPP_DE)
length(DE_genes)
F_notatus_BW = read.csv("F_notatus_DE_TR_FW.csv")
F_notatus_DE = unique(F_notatus_BW$GeneID)
length(F_notatus_DE)
DE_genes = union(DE_genes,F_notatus_DE)
length(DE_genes)
F_olivaceus_BW = read.csv("F_olivaceus_DE_TR_BW.csv")
F_olivaceus_DE = unique(F_olivaceus_BW$GeneID)
length(F_olivaceus_DE)
DE_genes = union(DE_genes,F_olivaceus_DE)
length(DE_genes)
F_parvapinis_BW = read.csv("F_parvapinis_DE_TR_BW.csv")
F_parvapinis_DE = unique(F_parvapinis_BW$GeneID)
length(F_parvapinis_DE)
DE_genes = union(DE_genes,F_parvapinis_DE)
length(DE_genes)
F_rathbuni_BW = read.csv("F_rathbuni_DE_TR_FW.csv")
F_rathbuni_DE = unique(F_rathbuni_BW$GeneID)
length(F_rathbuni_DE)
DE_genes = union(DE_genes,F_rathbuni_DE)
length(DE_genes)
F_similis_BW = read.csv("F_similis_DE_TR_BW.csv")
F_similis_DE = unique(F_similis_BW$GeneID)
length(F_similis_DE)
DE_genes = union(DE_genes,F_similis_DE)
length(DE_genes)
L_goodei_BW = read.csv("L_goodei_DE_TR_BW.csv")
L_goodei_DE = unique(L_goodei_BW$GeneID)
length(L_goodei_DE)
DE_genes = union(DE_genes,L_goodei_DE)
length(DE_genes)
L_parva_BW = read.csv("L_parva_DE_TR_BW.csv")
L_parva_DE = unique(L_parva_BW$GeneID)
length(L_parva_DE)
DE_genes = union(DE_genes,L_parva_DE)
length(DE_genes)
DE_genes_BW <- DE_genes

# FW,BW
sub_F_similis_FW<-F_similis_data[F_similis_data$GeneID %in% DE_genes_FW,]
sub_F_similis_FW<-sub_F_similis_FW[!duplicated(sub_F_similis_FW$GeneID),]
sub_F_similis_BW<-F_similis_data[F_similis_data$GeneID %in% DE_genes_BW,]
sub_F_similis_BW<-sub_F_similis_BW[!duplicated(sub_F_similis_BW$GeneID),]

sub_L_parva_FW<-L_parva_data[L_parva_data$GeneID %in% DE_genes_FW,]
sub_L_parva_FW<-sub_L_parva_FW[!duplicated(sub_L_parva_FW$GeneID),]
sub_L_parva_BW<-L_parva_data[L_parva_data$GeneID %in% DE_genes_BW,]
sub_L_parva_BW<-sub_L_parva_BW[!duplicated(sub_L_parva_BW$GeneID),]

sub_L_goodei_FW<-L_goodei_data[L_goodei_data$GeneID %in% DE_genes_FW,]
sub_L_goodei_FW<-sub_L_goodei_FW[!duplicated(sub_L_goodei_FW$GeneID),]
sub_L_goodei_BW<-L_goodei_data[L_goodei_data$GeneID %in% DE_genes_BW,]
sub_L_goodei_BW<-sub_L_goodei_BW[!duplicated(sub_L_goodei_BW$GeneID),]

sub_F_rathbuni_FW<-F_rathbuni_data[F_rathbuni_data$GeneID %in% DE_genes_FW,]
sub_F_rathbuni_FW<-sub_F_rathbuni_FW[!duplicated(sub_F_rathbuni_FW$GeneID),]
sub_F_rathbuni_BW<-F_rathbuni_data[F_rathbuni_data$GeneID %in% DE_genes_BW,]
sub_F_rathbuni_BW<-sub_F_rathbuni_BW[!duplicated(sub_F_rathbuni_BW$GeneID),]

sub_F_parvapinis_FW<-F_parvapinis_data[F_parvapinis_data$GeneID %in% DE_genes_FW,]
sub_F_parvapinis_FW<-sub_F_parvapinis_FW[!duplicated(sub_F_parvapinis_FW$GeneID),]
sub_F_parvapinis_BW<-F_parvapinis_data[F_parvapinis_data$GeneID %in% DE_genes_BW,]
sub_F_parvapinis_BW<-sub_F_parvapinis_BW[!duplicated(sub_F_parvapinis_BW$GeneID),]

sub_F_olivaceus_FW<-F_olivaceus_data[F_olivaceus_data$GeneID %in% DE_genes_FW,]
sub_F_olivaceus_FW<-sub_F_olivaceus_FW[!duplicated(sub_F_olivaceus_FW$GeneID),]
sub_F_olivaceus_BW<-F_olivaceus_data[F_olivaceus_data$GeneID %in% DE_genes_BW,]
sub_F_olivaceus_BW<-sub_F_olivaceus_BW[!duplicated(sub_F_olivaceus_BW$GeneID),]

sub_F_notatus_FW<-F_notatus_data[F_notatus_data$GeneID %in% DE_genes_FW,]
sub_F_notatus_FW<-sub_F_notatus_FW[!duplicated(sub_F_notatus_FW$GeneID),]
sub_F_notatus_BW<-F_notatus_data[F_notatus_data$GeneID %in% DE_genes_BW,]
sub_F_notatus_BW<-sub_F_notatus_BW[!duplicated(sub_F_notatus_BW$GeneID),]

sub_F_heteroclitusMDPP_FW<-F_heteroclitusMDPP_data[F_heteroclitusMDPP_data$GeneID %in% DE_genes_FW,]
sub_F_heteroclitusMDPP_FW<-sub_F_heteroclitusMDPP_FW[!duplicated(sub_F_heteroclitusMDPP_FW$GeneID),]
sub_F_heteroclitusMDPP_BW<-F_heteroclitusMDPP_data[F_heteroclitusMDPP_data$GeneID %in% DE_genes_BW,]
sub_F_heteroclitusMDPP_BW<-sub_F_heteroclitusMDPP_BW[!duplicated(sub_F_heteroclitusMDPP_BW$GeneID),]

sub_F_heteroclitusMDPL_FW<-F_heteroclitusMDPL_data[F_heteroclitusMDPL_data$GeneID %in% DE_genes_FW,]
sub_F_heteroclitusMDPL_FW<-sub_F_heteroclitusMDPL_FW[!duplicated(sub_F_heteroclitusMDPL_FW$GeneID),]
sub_F_heteroclitusMDPL_BW<-F_heteroclitusMDPL_data[F_heteroclitusMDPL_data$GeneID %in% DE_genes_BW,]
sub_F_heteroclitusMDPL_BW<-sub_F_heteroclitusMDPL_BW[!duplicated(sub_F_heteroclitusMDPL_BW$GeneID),]

sub_F_grandis_FW<-F_grandis_data[F_grandis_data$GeneID %in% DE_genes_FW,]
sub_F_grandis_FW<-sub_F_grandis_FW[!duplicated(sub_F_grandis_FW$GeneID),]
sub_F_grandis_BW<-F_grandis_data[F_grandis_data$GeneID %in% DE_genes_BW,]
sub_F_grandis_BW<-sub_F_grandis_BW[!duplicated(sub_F_grandis_BW$GeneID),]

sub_F_chrysotus_FW<-F_chrysotus_data[F_chrysotus_data$GeneID %in% DE_genes_FW,]
sub_F_chrysotus_FW<-sub_F_chrysotus_FW[!duplicated(sub_F_chrysotus_FW$GeneID),]
sub_F_chrysotus_BW<-F_chrysotus_data[F_chrysotus_data$GeneID %in% DE_genes_BW,]
sub_F_chrysotus_BW<-sub_F_chrysotus_BW[!duplicated(sub_F_chrysotus_BW$GeneID),]

sub_A_xenica_FW<-A_xenica_data[A_xenica_data$GeneID %in% DE_genes_FW,]
sub_A_xenica_FW<-sub_A_xenica_FW[!duplicated(sub_A_xenica_FW$GeneID),]
sub_A_xenica_BW<-A_xenica_data[A_xenica_data$GeneID %in% DE_genes_BW,]
sub_A_xenica_BW<-sub_A_xenica_BW[!duplicated(sub_A_xenica_BW$GeneID),]

# Combine all FW
all <- merge(sub_F_rathbuni_FW,sub_F_parvapinis_FW,by="GeneID")
all <- merge(all,sub_A_xenica_FW,by="GeneID")             
all <- merge(all,sub_F_chrysotus_FW,by="GeneID")             
all <- merge(all,sub_F_grandis_FW,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPL_FW,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPP_FW,by="GeneID")             
all <- merge(all,sub_F_notatus_FW,by="GeneID")             
all <- merge(all,sub_F_olivaceus_FW,by="GeneID")             
all <- merge(all,sub_L_goodei_FW,by="GeneID")             
all <- merge(all,sub_L_parva_FW,by="GeneID")
all <- merge(all,sub_F_similis_FW,by="GeneID")
# Separate gene lists out
all_DE_FW <- all
all_DE_FW<-all_DE_FW[,c(2:106)]
dim(all_DE_FW)
colnames(all_DE_FW)
d<-as.matrix(all_DE_FW)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="TR vs. FW DE genes", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
# Combine all BW
all <- merge(sub_F_rathbuni_BW,sub_F_parvapinis_BW,by="GeneID")
all <- merge(all,sub_A_xenica_BW,by="GeneID")             
all <- merge(all,sub_F_chrysotus_BW,by="GeneID")             
all <- merge(all,sub_F_grandis_BW,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPL_BW,by="GeneID")             
all <- merge(all,sub_F_heteroclitusMDPP_BW,by="GeneID")             
all <- merge(all,sub_F_notatus_BW,by="GeneID")             
all <- merge(all,sub_F_olivaceus_BW,by="GeneID")             
all <- merge(all,sub_L_goodei_BW,by="GeneID")             
all <- merge(all,sub_L_parva_BW,by="GeneID")
all <- merge(all,sub_F_similis_BW,by="GeneID")
# Separate gene lists out
all_DE_BW <- all
all_DE_BW<-all_DE_BW[,c(2:106)]
colnames(all_DE_BW)
d<-as.matrix(all_DE_BW)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="TR vs. BW DE genes", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
