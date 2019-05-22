
library(data.table)

setwd("~/Documents/UCDavis/Whitehead")
condition_DE <- read.csv("15_ppt_v_0.2_ppt.csv")
physiology_DE <-read.csv("M_v_FW.csv")
M <- read.csv("15_ppt_v_0.2_ppt_M.csv",stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
FW <- read.csv("15_ppt_v_0.2_ppt_FW.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade3_M <- read.csv("15_ppt_v_0.2_ppt_Clade3_M.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade3_FW <- read.csv("15_ppt_v_0.2_ppt_Clade3_FW.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade2_M <- read.csv("15_ppt_v_0.2_ppt_Clade2_M.csv",stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade2_FW <- read.csv("15_ppt_v_0.2_ppt_Clade2_FW.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade1_M <- read.csv("15_ppt_v_0.2_ppt_Clade1_M.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
Clade1_FW <- read.csv("15_ppt_v_0.2_ppt_Clade1_FW.csv", stringsAsFactors = FALSE, header = TRUE, row.names = NULL)

# sig genes 15ppt vs. 0.2ppt

dim(condition_DE)

condition_sig<-condition[condition$adj.P.Val <= 0.05,]
condition_sig<-condition_sig[!is.na(condition_sig$adj.P.Val),]
condition_sig<-condition_sig$X
length(condition_sig)

# sig genes M vs. FW
dim(physiology_DE)

physiology_sig<-physiology_DE[physiology_DE$adj.P.Val <= 0.05,]
physiology_sig<-physiology_sig[!is.na(physiology_sig$adj.P.Val),]
physiology_sig<-physiology_sig$X
length(physiology_sig)

# sig genes M and FW
# two-way
# conserved response
dim(M)
M_sig <- M[M$adj.P.Val <= 0.05,] 
M_sig <- M_sig[!is.na(M_sig$adj.P.Val),]
M_sig <- M_sig$X
length(M_sig)
dim(FW)
FW_sig <- FW[FW$adj.P.Val <= 0.05,] 
FW_sig <- FW_sig[!is.na(FW_sig$adj.P.Val),]
FW_sig <- FW_sig$X
length(FW_sig)

# Clade 3 - specific response
# sig genes Clade 3, M
# but not sig for 2-way

Clade3_M_sig <- Clade3_M[Clade3_M$adj.P.Val <= 0.05,]
Clade3_M_sig <- Clade3_M_sig[!is.na(Clade3_M_sig$adj.P.Val),]
dim(Clade3_M_sig)
Clade3_M_sig_specific <- Clade3_M_sig[!Clade3_M_sig$X %in% M_sig,]
dim(Clade3_M_sig_specific)
Clade3_M_sig_specific <- Clade3_M_sig_specific$X

# Clade 3 - specific response
# sig genes Clade 1,FW
# but not sig for 2-way

Clade3_FW_sig <- Clade3_FW[Clade3_FW$adj.P.Val <= 0.05,]
Clade3_FW_sig <- Clade3_FW_sig[!is.na(Clade3_FW_sig$adj.P.Val),]
dim(Clade3_FW_sig)
Clade3_FW_sig_specific <- Clade3_FW_sig[!Clade3_FW_sig$X %in% FW_sig,]
dim(Clade3_FW_sig_specific)
Clade3_FW_sig_specific <- Clade3_FW_sig_specific$X

# Clade 2 - specific response
# sig genes Clade 2, M
# but not sig for 2-way

Clade2_M_sig <- Clade2_M[Clade2_M$adj.P.Val <= 0.05,]
Clade2_M_sig <- Clade2_M_sig[!is.na(Clade2_M_sig$adj.P.Val),]
dim(Clade2_M_sig)
Clade2_M_sig_specific <- Clade2_M_sig[!Clade2_M_sig$X %in% M_sig,]
dim(Clade2_M_sig_specific)
Clade2_M_sig_specific <- Clade2_M_sig_specific$X

# Clade 2 - specific response
# sig genes Clade 2, FW
# but not sig for 2-way

Clade2_FW_sig <- Clade2_FW[Clade2_FW$adj.P.Val <= 0.05,]
Clade2_FW_sig <- Clade2_FW_sig[!is.na(Clade2_FW_sig$adj.P.Val),]
dim(Clade2_FW_sig)
Clade2_FW_sig_specific <- Clade2_FW_sig[!Clade2_FW_sig$X %in% FW_sig,]
dim(Clade2_FW_sig_specific)
Clade2_FW_sig_specific <- Clade2_FW_sig_specific$X

# Clade 1 - specific response
# sig genes Clade 1, M
# but not sig for 2-way

Clade1_M_sig <- Clade1_M[Clade1_M$adj.P.Val <= 0.05,]
Clade1_M_sig <- Clade1_M_sig[!is.na(Clade1_M_sig$adj.P.Val),]
dim(Clade1_M_sig)
Clade1_M_sig_specific <- Clade1_M_sig[!Clade1_M_sig$X %in% M_sig,]
dim(Clade1_M_sig_specific)
Clade1_M_sig_specific <- Clade1_M_sig_specific$X

# Clade 1 - specific response
# sig genes Clade 1, FW
# but not sig for 2-way

Clade1_FW_sig <- Clade1_FW[Clade1_FW$adj.P.Val <= 0.05,]
Clade1_FW_sig <- Clade1_FW_sig[!is.na(Clade1_FW_sig$adj.P.Val),]
dim(Clade1_FW_sig)
Clade1_FW_sig_specific <- Clade1_FW_sig[!Clade1_FW_sig$X %in% FW_sig,]
dim(Clade1_FW_sig_specific)
Clade1_FW_sig_specific <- Clade1_FW_sig_specific$X

# if sig for two-way, then conserved response 
# but not clade-specific


# take averages of normalized expression for each treatment, for each species, 

norm_counts <- read.csv("normalized_counts.csv",stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
colnames(norm_counts)
cols.norm_counts <- colnames(norm_counts)

# F_rathbuni - Clade 1 - FW
# F_grandis - Clade 1 - M
# F_notatus - Clade 3 - FW
# F_parvapinis - Clade 2 - M
# L_goodei - Clade 2 - FW
# F_olivaceus - Clade 3 - FW
# L_parva - Clade 2 - M
# F_heteroclitusMDPP - Clade 1 - M
# F_similis - Clade 1 - M
# F_diaphanus - Clade 1 - M
# F_chrysotus - Clade 3 - M
# A_xenica - Clade 3 - M
# F_catanatus - Clade 1 - FW
# F_heteroclitusMDPL - Clade 1 - M

species_condition<-as.vector(paste(species,condition,sep="."))
rownames(norm_counts) <- norm_counts$Ensembl
norm_counts <- norm_counts[,-1]
colnames(norm_counts) <- species_condition
dim(norm_counts) 
mean_norm_counts<-sapply(unique(colnames(norm_counts)), function(i)
  rowMeans(norm_counts[,colnames(norm_counts) == i]))

# F_grandis - Clade 1 - M
# F_similis - Clade 1 - M
# F_diaphanus - Clade 1 - M
# F_heteroclitusMDPL - Clade 1 - M
# F_heteroclitusMDPP - Clade 1 - M
# F_catanatus - Clade 1 - FW
# F_rathbuni - Clade 1 - FW
# F_parvapinis - Clade 2 - M
# L_parva - Clade 2 - M
# L_goodei - Clade 2 - FW
# F_chrysotus - Clade 3 - M
# A_xenica - Clade 3 - M
# F_notatus - Clade 3 - FW
# F_olivaceus - Clade 3 - FW

ph <- c("M","M","M","M","M","M","M","M","M","M","FW","FW","FW","FW","M","M","M","M","FW","FW","M","M","M","M","FW","FW","FW","FW")
cl <- c("Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade2","Clade2","Clade2","Clade2","Clade2","Clade2","Clade3","Clade3","Clade3","Clade3","Clade3","Clade3","Clade3","Clade3")
condition <- c("15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt")
  
colnames(mean_norm_counts)
sample_order <- c("F_grandis.15_ppt","F_grandis.0.2_ppt","F_similis.15_ppt","F_similis.0.2_ppt",
           "F_diaphanus.15_ppt","F_diaphanus.0.2_ppt","F_heteroclitusMDPL.15_ppt","F_heteroclitusMDPL.0.2_ppt",
           "F_heteroclitusMDPP.15_ppt","F_heteroclitusMDPP.0.2_ppt","F_catanatus.15_ppt","F_catanatus.0.2_ppt",
           "F_rathbuni.15_ppt","F_rathbuni.0.2_ppt","F_parvapinis.15_ppt","F_parvapinis.0.2_ppt",
           "L_parva.15_ppt","L_parva.0.2_ppt","L_goodei.15_ppt","L_goodei.0.2_ppt",
           "F_chrysotus.15_ppt","F_chrysotus.0.2_ppt","A_xenica.15_ppt","A_xenica.0.2_ppt",
           "F_notatus.15_ppt","F_notatus.0.2_ppt","F_olivaceous.15_ppt","F_olivaceous.0.2_ppt")

mean_norm_counts_ordered <- mean_norm_counts[,sample_order]
colnames(mean_norm_counts_ordered)

length(Clade1_M_sig_specific)
mean_norm_counts_ordered_M_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% Clade1_M_sig_specific,]
dim(mean_norm_counts_ordered_M_sig)

length(FW_sig)
mean_norm_counts_ordered_FW_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% FW_sig,]

d<-as.matrix(mean_norm_counts_ordered_M_sig)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="M, Clade1-specific, padj<0.05)",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)

mean_norm_counts_ordered_condition_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% condition_sig,]
dim(mean_norm_counts_ordered_condition_sig)

rld <- log2(mean_norm_counts_ordered_condition_sig+1)
geneDists <- dist(mean_norm_counts_ordered_condition_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = geneDists, 
         cluster_cols= FALSE,
         annotation_col=df,
         scale = "row")


mean_norm_counts_ordered_M_sig <- mean_norm_counts_ordered_M_sig[rownames(mean_norm_counts_ordered_M_sig) %in% condition_sig,]

mean_norm_counts_ordered_physiology_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% physiology_sig,]
dim(mean_norm_counts_ordered_physiology_sig)
rld <- log2(mean_norm_counts_ordered_physiology_sig+1)
geneDists <- dist(mean_norm_counts_ordered_physiology_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = geneDists, cluster_cols= FALSE,annotation_col=df)
