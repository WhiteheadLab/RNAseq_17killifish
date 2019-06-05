library(gridExtra)
library(data.table)
library(dplyr)
library(grid)
library(ggplot2)
library(lattice)
library(biomaRt)
library(RColorBrewer)
library(viridis)

setwd("~/Documents/UCDavis/Whitehead")
# normalized counts
norm_counts <- read.csv("normalized_counts.csv",stringsAsFactors = FALSE, header = TRUE,row.names = "Ensembl")


# ============================================
# biomart annotation
# https://uswest.ensembl.org/Fundulus_heteroclitus/Info/Index
# ============================================

ensembl=useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("fheteroclitus_gene_ensembl",mart=ensembl)
ensembl_proteinID = rownames(norm_counts)
length(ensembl_proteinID)
ann<-getBM(attributes=c('ensembl_peptide_id','ensembl_transcript_id','ensembl_gene_id','gene_biotype','external_gene_name','description','entrezgene'), filters = 'ensembl_peptide_id', values = ensembl_proteinID, mart=ensembl)
head(ann)
dim(ann)
length(unique(ann$ensembl_peptide_id))
ann <- ann[!duplicated(ann[,c(1)]),]
dim(ann)

#-------------------------------------
# contrasts
#-------------------------------------

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

#-------------------------------------
# main and interaction effects, F tests
#-------------------------------------

main_condition <- read.csv("main_condition.csv", stringsAsFactors = FALSE,header = TRUE)

main_phys <- read.csv("main_phys.csv", stringsAsFactors = FALSE, header = TRUE)
main_clade <- read.csv("main_clade.csv", stringsAsFactors = FALSE, header = TRUE)


int_phys_condition <- read.csv("interaction_physiology_condition.csv", stringsAsFactors = FALSE, header = TRUE)
int_phys_clade <- read.csv("interaction_physiology_clade.csv", stringsAsFactors = FALSE, header = TRUE)
int_clade_condition <- read.csv("interaction_clade_condition.csv", stringsAsFactors = FALSE, header = TRUE)
int_three <- read.csv("interaction_threeway.csv", stringsAsFactors = FALSE, header = TRUE)

#-------------------------------------
# contrasts
# sig genes 15ppt vs. 0.2ppt
#-------------------------------------

dim(condition_DE)
condition_sig<-condition[condition$adj.P.Val <= 0.05,]
condition_sig<-condition_sig[!is.na(condition_sig$adj.P.Val),]
condition_sig<-condition_sig$X
length(condition_sig)

tmp <- norm_counts[rownames(norm_counts) %in% condition_sig,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
rownames(tmp_ann) <- tmp_ann$Row.names
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()


pdf("~/Documents/UCDavis/Whitehead/15ppt_v_0.2ppt_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

# sig genes M vs. FW
physiology_sig<-physiology_DE[physiology_DE$adj.P.Val <= 0.05,]
physiology_sig<-physiology_sig[!is.na(physiology_sig$adj.P.Val),]
physiology_sig<-physiology_sig$X
length(physiology_sig)

tmp <- norm_counts[rownames(norm_counts) %in% physiology_sig,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
rownames(tmp_ann) <- tmp_ann$Row.names
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()


pdf("~/Documents/UCDavis/Whitehead/physiology_M_v_FW_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

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
Clade3_M_sig_genes <- Clade3_M_sig$X
length(Clade3_M_sig_genes)
[1] 105

Clade3_M_sig_specific <- Clade3_M_sig[!Clade3_M_sig$X %in% M_sig,]
dim(Clade3_M_sig_specific)
Clade3_M_sig_specific <- Clade3_M_sig_specific$X
length(Clade3_M_sig_specific)

length(Clade3_M_sig_specific)
# [1] 50
tmp <- norm_counts[rownames(norm_counts) %in% Clade3_M_sig_specific,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann <- tmp_ann[!(tmp_ann$external_gene_name == ""),]
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()


pdf("~/Documents/UCDavis/Whitehead/clade3_M_specific_conditioncladephysiology_15ppt_v_0.2ppt_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

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

#-------------------------------------
# significant genes
# main effects
#-------------------------------------

sig_main_phys <- main_phys[main_phys$adj.P.Val<0.05,]
sig_main_phys <- sig_main_phys$X
length(sig_main_phys)
sig_main_condition <- main_condition[main_condition$adj.P.Val<0.05,]
sig_main_condition <- sig_main_condition$X
length(sig_main_condition)
sig_main_clade <- main_clade[main_clade$adj.P.Val<0.05,]
sig_main_clade <- sig_main_clade$X
length(sig_main_clade)
# two-way interactions
sig_int_phys_condition <- int_phys_condition[int_phys_condition$adj.P.Val<0.05,]
sig_int_phys_condition <- sig_int_phys_condition$X
length(sig_int_phys_condition)
sig_int_phys_clade <- int_phys_clade[int_phys_clade$adj.P.Val<0.05,]
sig_int_phys_clade <- sig_int_phys_clade$X
length(sig_int_phys_clade)
sig_int_clade_condition <- int_clade_condition[int_clade_condition$adj.P.Val<0.05,]
sig_int_clade_condition <- sig_int_clade_condition$X
length(sig_int_clade_condition)
# three-way interaction
sig_int_three <- int_three[int_three$adj.P.Val<0.05,]
sig_int_three <- sig_int_three$X
length(sig_int_three)

# --------------------------
# import normalized counts
# --------------------------
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
condition <- c("0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt","0.2_ppt","15_ppt")
  
colnames(mean_norm_counts)
sample_order <- c("F_grandis.0.2_ppt","F_grandis.15_ppt","F_similis.0.2_ppt","F_similis.15_ppt","F_diaphanus.0.2_ppt",
           "F_diaphanus.15_ppt","F_heteroclitusMDPL.0.2_ppt","F_heteroclitusMDPL.15_ppt",
           "F_heteroclitusMDPP.0.2_ppt","F_heteroclitusMDPP.15_ppt","F_catanatus.0.2_ppt","F_catanatus.15_ppt",
           "F_rathbuni.0.2_ppt","F_rathbuni.15_ppt","F_parvapinis.0.2_ppt","F_parvapinis.15_ppt",
           "L_parva.0.2_ppt","L_parva.15_ppt","L_goodei.0.2_ppt","L_goodei.15_ppt",
           "F_chrysotus.0.2_ppt","F_chrysotus.15_ppt","A_xenica.0.2_ppt","A_xenica.15_ppt",
           "F_notatus.0.2_ppt","F_notatus.15_ppt","F_olivaceous.0.2_ppt","F_olivaceous.15_ppt")

mean_norm_counts_ordered <- mean_norm_counts[,sample_order]
colnames(mean_norm_counts_ordered)

length(sig_main_phys)
length(sig_main_condition)
length(sig_main_clade)
# two-way interactions
length(sig_int_phys_condition)
length(sig_int_phys_clade)
length(sig_int_clade_condition)
# three-way interaction
length(sig_int_three)

mean_norm_counts_ordered_phys_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_phys,]
dim(mean_norm_counts_ordered_phys_sig)
mean_norm_counts_ordered_cond_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_condition,]
dim(mean_norm_counts_ordered_cond_sig)
mean_norm_counts_ordered_clade_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_clade,]
dim(mean_norm_counts_ordered_clade_sig)


d<-as.matrix(mean_norm_counts_ordered_phys_sig)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="condition main effects, padj<0.05)",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)

rld <- log2(mean_norm_counts_ordered_phys_sig+1)
geneDists <- dist(mean_norm_counts_ordered_phys_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = geneDists, 
         cluster_cols= FALSE,
         annotation_col=df,
         scale = "row")
head(rld)

# for each FW col, subtract values by itself and subtract 15ppt col by original FW
# in for loop

# First, try this on two cols
# make new df
colnames(rld)
odd <- seq(1,ncol(rld),2)
even <- seq(2,ncol(rld),2)
rld_sub <- data.frame(matrix("", ncol = length(colnames(rld)), nrow = length(rownames(rld))))  
for (i in 1:14){
  print(i)
  print(odd[i])
  print(even[i])
  rld_sub <- cbind(rld_sub,rld[,odd[i]] - rld[,odd[i]])
  rld_sub <- cbind(rld_sub,rld[,even[i]] - rld[,odd[i]])
}
rld_sub <- rld_sub[,c(29:56)]
colnames(rld_sub) <- colnames(rld)
dim(rld_sub)
geneDists <- dist(rld_sub)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld_sub)
pheatmap(rld_sub, show_rownames=FALSE,
         clustering_distance_rows = geneDists, 
         cluster_cols= FALSE,
         annotation_col=df,
         scale = "none",
         color = colorRampPalette(c("mediumblue", "white", "yellow"))(20))

head(rld)







mean_norm_counts_ordered_M_sig <- mean_norm_counts_ordered_M_sig[rownames(mean_norm_counts_ordered_M_sig) %in% condition_sig,]

mean_norm_counts_ordered_physiology_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% physiology_sig,]
dim(mean_norm_counts_ordered_physiology_sig)
rld <- log2(mean_norm_counts_ordered_physiology_sig+1)
geneDists <- dist(mean_norm_counts_ordered_physiology_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = geneDists, cluster_cols= FALSE,annotation_col=df,scale="row")

# ----------------------

# plotting expression
# of GOI

all_goi<-c("ENSFHEP00000003519","ENSFHEP00000011077","ENSFHEP00000015386","ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
           "ENSFHEP00000015383","ENSFHEP00000009753","ENSFHEP00000006725","ENSFHEP00000008393",
           "ENSFHEP00000013324","ENSFHEP00000001609","ENSFHEP00000013324","ENSFHEP00000034177",
           "ENSFHEP00000015765","ENSFHEP00000017303","ENSFHEP00000000036","ENSFHEP00000031108",
           "ENSFHEP00000016853","ENSFHEP00000003908","ENSFHEP00000026411","ENSFHEP00000006282")
tmp <- norm_counts[rownames(norm_counts) %in% all_goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)

tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene) +
  labs(x="Clade:Condition",
       y="Expression (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1))

# ----------------------------
# main effects
# Physiology
# ----------------------------
length(sig_main_phys_genes)

tmp <- norm_counts[rownames(norm_counts) %in% sig_main_phys_genes,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$Row.names

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
          aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
        }
dev.off()




# ----------------------------
# main effects
# Condition
# ----------------------------
length(sig_main_condition_genes)
tmp <- norm_counts[rownames(norm_counts) %in% sig_main_condition_genes,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann <- tmp_ann[!(tmp_ann$external_gene_name == ""),]
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene) +
  labs(x="Clade:Condition",
       y="Expression (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1))


# ----------------------------
# main effects
# Clade
# ----------------------------

length(sig_main_clade_genes)
# [1] 2614
tmp <- j[rownames(j) %in% sig_main_clade_genes,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/clade_maineffect_28May2019_2.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()


# ---------------------------
# two-way interactions
length(sig_int_phys_condition_genes)
tmp <- norm_counts[rownames(norm_counts) %in% sig_int_phys_condition_genes,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann <- tmp_ann[!(tmp_ann$external_gene_name == ""),]
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/physiologycondition_2wayinteraction_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

length(sig_int_phys_clade_genes)
tmp <- j[rownames(j) %in% sig_int_phys_clade_genes,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/physiologyclade_2wayinteraction_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()


length(sig_int_clade_condition_genes)

tmp <- norm_counts[rownames(norm_counts) %in% sig_int_clade_condition_genes,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann <- tmp_ann[!(tmp_ann$external_gene_name == ""),]
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/cladecondition_2wayinteraction_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

# --------------------------
# three-way interaction
length(sig_int_three_genes)
sig_int_three_genes

tmp <- norm_counts[rownames(norm_counts) %in% sig_int_three_genes,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
rownames(tmp_ann) <- tmp_ann$external_gene_name
tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

pdf("~/Documents/UCDavis/Whitehead/threewayinteraction_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 20)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
dev.off()

# ----------------------------
# MA plot


gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")

gene_id <- c("avpr2aa","slc24a5","CLDN4","aqp3","cftr","kcnj2a","polyamine-modulated factor 1-like","kcnj1a.6","sept2B","septin-2", "cipcb","clcn2c","zymogen granule membrane protein 16","atp1a1b","solute carrier family 24 member 2")
protein_id <- c("ENSFHEP00000000036","ENSFHEP00000001609","ENSFHEP00000003908","ENSFHEP00000006725","ENSFHEP00000008393","ENSFHEP00000009753","ENSFHEP00000013324",
                "ENSFHEP00000015383","ENSFHEP00000015765","ENSFHEP00000016853","ENSFHEP00000017303","ENSFHEP00000019510","ENSFHEP00000025841",
                "ENSFHEP00000031108","ENSFHEP00000034177")
res<-condition_DE
res_ordered <-as.data.frame(condition_DE[order(condition_DE$adj.P.Val),])
res_filtered <-subset(res_ordered,res_ordered$padj<0.05)
id<-res_filtered$X
res_filtered<-cbind(res_filtered,id)
plot(res$AveExpr, res$logFC, 
     col=ifelse(res$adj.P.Val < 0.05, "red","gray67"),
     main="15_ppt vs. 0.2_ppt (padj<0.05)",pch=20,cex=1)
abline(h=c(-1,1), col="blue")

resSig = res_ordered[res_ordered$X %in% protein_id, ]
dim(resSig)
genes<-rownames(resSig)
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"AveExpr"]
log2FoldChange_mygenes <- mygenes[,"logFC"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=gene_id,pos=2,cex=0.60)



# ----------------
# single genes plots
# evolutionary pattern examples, genes of interest
# ------------------


# ----------------------------
# main effects
# Physiology
# ----------------------------
length(sig_main_phys_genes)
# CFTR
goi <- c("ENSFHEP00000008393")

tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# Na/K ATPases
goi<-c("ENSFHEP00000003519","ENSFHEP00000011077","ENSFHEP00000015386","ENSFHEP00000031108")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000001223
# DOT1 like histone lysine methyltransferase
goi<-c("ENSFHEP00000001223")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()


# ENSFHEP00000006401
# This looks like complete crap.
# for some reason this has sig effect of condition
# low expressed - therefore picked up as significant
goi<-c("ENSFHEP00000006401")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000014138
# for some reason this has a sig effect of condition
# myosin heavy chain
goi<-c("ENSFHEP00000014138")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$Row.names

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000020413
# claudin 10
goi<-c("ENSFHEP00000020413")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# sig for three-way
# ENSFHEP00000023795
# thyroid hormone
goi<-c("ENSFHEP00000023795")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000027492
goi<-c("ENSFHEP00000027492")
# gfi1ab
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()


# ENSFHEP00000031409
goi<-c("ENSFHEP00000031409")
# gfi1ab
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000011721
goi<-c("ENSFHEP00000011721")
# sterol 26-hydroxylase
# cyp27a7
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()

# ENSFHEP00000033747
# per3
goi<-c("ENSFHEP00000033747")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}

# ENSFHEP00000023833
# nocta 
# nocturnin
goi<-c("ENSFHEP00000023833")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}



goi<-c("ENSFHEP00000000225")
tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann) <- tmp_ann$external_gene_name

tmp <- tmp_ann[,c(2:82)]
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_28May2019.pdf",paper="USr",width=13.5, height=8)
for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ], 
               aes(x=clade:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)",
               fill = "Native Physiology") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}

# --------------------
# faceted Log2FC correlations
# for each clade, 
# M, 15 vs 0.2ppt
# FW, 15 vs. 0.2ppt
# --------------------

# Clade 3, M
condition_sig
Clade3_M_sig<-Clade3_M[Clade3_M$X %in% condition_sig,]
Clade3_M_sig<-Clade3_M_sig[!is.na(Clade3_M_sig$adj.P.Val),]
Clade3_M_sig_genes<-Clade3_M_sig$X
length(Clade3_M_sig_genes)
clade <- data.frame(factor(rep("Clade3",length(Clade3_M_sig_genes))))
physiology <- data.frame(factor(rep("M",length(Clade3_M_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade3_M_sig$X
tmp <- data.frame(Clade3_M_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade3_M_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade3_M_logFC <- tcounts
dim(Clade3_M_logFC)

# Clade 3, FW
Clade3_FW_sig<-Clade3_FW[Clade3_FW$X %in% condition_sig,]
Clade3_FW_sig<-Clade3_FW_sig[!is.na(Clade3_FW_sig$adj.P.Val),]
Clade3_FW_sig_genes<-Clade3_FW_sig$X
length(Clade3_FW_sig_genes)
clade <- data.frame(factor(rep("Clade3",length(Clade3_FW_sig_genes))))
physiology <- data.frame(factor(rep("FW",length(Clade3_FW_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade3_FW_sig$X
tmp <- data.frame(Clade3_FW_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade3_FW_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade3_FW_logFC <- tcounts
dim(Clade3_FW_logFC)
dim(Clade3_M_logFC)
Clade_Phys_logFC <- rbind(Clade3_M_logFC,Clade3_FW_logFC)
dim(Clade_Phys_logFC)
x<- spread(Clade_Phys_logFC, Physiology, logFC)

# Clade 2, M
Clade2_M_sig<-Clade2_M[Clade2_M$X %in% condition_sig,,]
Clade2_M_sig<-Clade2_M_sig[!is.na(Clade2_M_sig$adj.P.Val),]
Clade2_M_sig_genes<-Clade2_M_sig$X
length(Clade2_M_sig_genes)
clade <- data.frame(factor(rep("Clade2",length(Clade2_M_sig_genes))))
physiology <- data.frame(factor(rep("M",length(Clade2_M_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade2_M_sig$X
tmp <- data.frame(Clade2_M_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade2_M_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade2_M_logFC <- tcounts
dim(Clade2_M_logFC)
Clade_Phys_logFC <- rbind(Clade_Phys_logFC,Clade2_M_logFC)
dim(Clade_Phys_logFC)
x<- spread(Clade_Phys_logFC, Physiology, logFC)
# Clade 2, FW
Clade2_FW_sig<-Clade2_FW[Clade2_FW$X %in% condition_sig,,]
Clade2_FW_sig<-Clade2_FW_sig[!is.na(Clade2_FW_sig$adj.P.Val),]
Clade2_FW_sig_genes<-Clade2_FW_sig$X
length(Clade2_FW_sig_genes)
clade <- data.frame(factor(rep("Clade2",length(Clade2_FW_sig_genes))))
physiology <- data.frame(factor(rep("FW",length(Clade2_FW_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade2_FW_sig$X
tmp <- data.frame(Clade2_FW_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade2_FW_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade2_FW_logFC <- tcounts
dim(Clade2_FW_logFC)
Clade_Phys_logFC <- rbind(Clade_Phys_logFC,Clade2_FW_logFC)
dim(Clade_Phys_logFC)
x<- spread(Clade_Phys_logFC, Physiology, logFC)

# Clade 1, M
Clade1_M_sig<-Clade1_M[Clade1_M$X %in% condition_sig,]
Clade1_M_sig<-Clade1_M_sig[!is.na(Clade1_M_sig$adj.P.Val),]
Clade1_M_sig_genes<-Clade1_M_sig$X
length(Clade1_M_sig_genes)
clade <- data.frame(factor(rep("Clade1",length(Clade1_M_sig_genes))))
physiology <- data.frame(factor(rep("M",length(Clade1_M_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade1_M_sig$X
tmp <- data.frame(Clade1_M_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade1_M_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade1_M_logFC <- tcounts
dim(Clade1_M_logFC)
Clade_Phys_logFC <- rbind(Clade_Phys_logFC,Clade1_M_logFC)
dim(Clade_Phys_logFC)
x<- spread(Clade_Phys_logFC, Physiology, logFC)

# Clade 1, FW
Clade1_FW_sig<-Clade1_FW[Clade1_FW$X %in% condition_sig,]
Clade1_FW_sig<-Clade1_FW_sig[!is.na(Clade1_FW_sig$adj.P.Val),]
Clade1_FW_sig_genes<-Clade1_FW_sig$X
length(Clade1_FW_sig_genes)
clade <- data.frame(factor(rep("Clade1",length(Clade1_FW_sig_genes))))
physiology <- data.frame(factor(rep("FW",length(Clade1_FW_sig_genes))))
info <- cbind(clade,physiology)
colnames(info) <- c("Clade","Physiology")
rownames(info) <- Clade1_FW_sig$X
tmp <- data.frame(Clade1_FW_sig$logFC)
colnames(tmp) <- c("logFC")
rownames(tmp) <- Clade1_FW_sig$X
tmp <- as.matrix(tmp)
tcounts <- tmp %>%
  merge(info, ., by="row.names")
Clade1_FW_logFC <- tcounts
dim(Clade1_FW_logFC)
Clade_Phys_logFC <- rbind(Clade_Phys_logFC,Clade1_FW_logFC)
dim(Clade_Phys_logFC)
colnames(Clade_Phys_logFC)
x <- spread(Clade_Phys_logFC, Physiology, logFC)

ggplot(x, aes(x = FW, y = M)) +
  geom_point(alpha = 0.5, size = 3) +
  facet_wrap(~clade_f,scales = "free") +
  theme_bw() +
  scale_color_viridis(name = "mean LFC",
                      breaks = c(-9, -6, -3, 0, 3,  6, 9),
                      limits=c(-9, 9))  +
  labs(#subtitle = "Comparison of log fold changes (LFC)",
    y = "LFC of genes from native M physiology, sig response 15 ppt vs. 0.2 ppt",
    x = "LFC of genes from native FW physiology, sig response 15 ppt vs. 0.2 ppt")

  