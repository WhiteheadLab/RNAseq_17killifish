library(gridExtra)
library(data.table)
library(dplyr)
library(grid)
library(ggplot2)
library(lattice)
library(biomaRt)
library(RColorBrewer)
library(viridis)
library(pheatmap)

setwd("~/Documents/UCDavis/Whitehead/kfish_expression_July2019/")
# normalized counts
norm_counts <- read.csv("normalized_counts.csv",stringsAsFactors = FALSE, header = TRUE,row.names = "Ensembl")
dim(norm_counts)

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
# significant genes
# main effects
#-------------------------------------
main_phys <- main_physiology
sig_main_phys <- main_phys[main_phys$adj.P.Val<0.05,]
sig_main_phys <- rownames(sig_main_phys)
length(sig_main_phys)
sig_main_condition <- main_salinity[main_salinity$adj.P.Val<0.05,]
sig_main_condition <- rownames(sig_main_salinity)
length(sig_main_condition)
sig_main_clade <- main_clade[main_clade$adj.P.Val<0.05,]
sig_main_clade <- rownames(sig_main_clade)
length(sig_main_clade)

# two-way interactions
dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)

int_phys_condition<- sig_salinity_physiology_interaction
sig_int_phys_condition <- int_phys_condition[int_phys_condition$adj.P.Val<0.05,]
sig_int_phys_condition <- rownames(sig_int_phys_condition)
length(sig_int_phys_condition)

int_phys_clade<- sig_clade_physiology_interaction
sig_int_phys_clade <- int_phys_clade[int_phys_clade$adj.P.Val<0.05,]
sig_int_phys_clade <- rownames(sig_int_phys_clade)
length(sig_int_phys_clade)

int_clade_condition<- sig_salinity_clade_interaction
sig_int_clade_condition <- int_clade_condition[int_clade_condition$adj.P.Val<0.05,]
sig_int_clade_condition <- rownames(sig_int_clade_condition)
length(sig_int_clade_condition)


# three-way interaction
int_three <- sig_threeway
sig_int_three <- int_three[int_three$adj.P.Val<0.05,]
sig_int_three <- rownames(sig_int_three)
length(sig_int_three)

# --------------------------
# import normalized counts
# formatting
# --------------------------
#norm_counts <- read.csv("normalized_counts.csv",stringsAsFactors = FALSE, header = TRUE, row.names = NULL)
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



# -----------------------------------

#-------------------------------------
# significant genes
#-------------------------------------

dim(sig_clade_physiology_interaction)
dim(sig_salinity_physiology_interaction)
dim(sig_salinity_clade_interaction)
dim(sig_threeway)
dim(sig_main_clade)
dim(sig_main_physiology)
dim(sig_main_salinity)


# sig expression
mean_norm_counts_ordered_phys_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_phys,]
dim(mean_norm_counts_ordered_phys_sig)
mean_norm_counts_ordered_cond_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_condition,]
dim(mean_norm_counts_ordered_cond_sig)
mean_norm_counts_ordered_clade_sig <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% sig_main_clade,]
dim(mean_norm_counts_ordered_clade_sig)

threeway <- rownames(sig_threeway)
main_clade <- rownames(sig_main_clade)
main_phys <- rownames(sig_main_physiology)
main_salinity <- rownames(sig_main_salinity)
int_phys_clade <- rownames(sig_clade_physiology_interaction)
int_salinity_phys <- rownames(sig_salinity_physiology_interaction)
int_salinity_clade <- rownames(sig_salinity_clade_interaction)

mean_threeway <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% threeway,]
mean_main_clade <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% main_clade,]
mean_main_phys <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% main_phys,]
mean_salinity <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% main_salinity,]
mean_phys_clade <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% int_phys_clade,]
mean_salinity_phys <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% int_salinity_phys,]
mean_salinity_clade <- mean_norm_counts_ordered[rownames(mean_norm_counts_ordered) %in% int_salinity_clade,]

mean_salinity_notinteraction_phys_sal <- mean_salinity[!rownames(mean_salinity) %in% int_salinity_phys,]

d<-as.matrix(mean_norm_counts_ordered_cond_sig)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="Interaction - physiology x salinity, padj<0.05)",
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)


rld <- log2(mean_norm_counts_ordered_cond_sig)
geneDists <- dist(mean_norm_counts_ordered_cond_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = "euclidean", 
         cluster_cols= FALSE,
         annotation_col=df,
         scale = "row")

# -------------------
# only BW samples
# -------------------
colnames(mean_norm_counts_ordered_cond_sig)
salinity_BW <- mean_norm_counts_ordered_cond_sig[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
rld <- log2(salinity_BW+1)
geneDists <- dist(salinity_BW)
ph_BW <- c("M","M","M","M","M","FW","FW","M","M","FW","M","M","FW","FW")
cl_BW <- c("Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade1","Clade2","Clade2","Clade2","Clade3","Clade3","Clade3","Clade3")
condition_BW <- c("15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt","15_ppt")
df <- data.frame(ph_BW,cl_BW, condition_BW,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = "euclidean", 
         cluster_cols= TRUE,
         annotation_col=df,
         scale = "row")
head(rld)


d<-rld
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



# for each FW col, subtract values by itself and subtract 15ppt col by original FW
# in for loop

# First, two cols
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
         clustering_distance_rows = "euclidean", 
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
tcounts %>% dplyr::select(Row.names, clade, species,physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

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
tcounts %>% dplyr::select(Row.names, clade, species, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

for (i in seq(1, length(unique(tcounts$gene)), 100)) {
  print(ggplot(tcounts[tcounts$gene %in% levels(as.factor(tcounts$gene))[i:(i+19)], ],
               aes(x=species:condition, y=expression,fill=physiology)) + 
          geom_boxplot() +
          facet_wrap(~gene) +
          labs(x="Clade:Condition",
               y="Expression (log2 cpm normalized counts)") +
          theme_bw() +
          theme(axis.text.x=element_text(angle=90, hjust=1)))
}
#dev.off()


C1<-ggplot(tcounts %>% filter(clade=='Clade1'), aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(size = 0.7) +
  ylim(-2.4, 10) +
  stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
  facet_grid(~species,scales='fixed', space = "free") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
  stat_summary(fun.y="mean", geom="point", size=5, 
               aes(shape = physiology,color=physiology)) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size= 20),
          axis.title.y = element_text(size= 20),
          axis.text.x = element_text(angle = 90, hjust = 1,size=20),
          strip.text.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 12)) +
    labs(y="Expression (log2 normalized counts)") +
    ggtitle("Clade 1")
plot(C1)
C2<-ggplot(tcounts %>% filter(clade=='Clade2'), aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
  ylim(-2.4, 10) +
    geom_point(size = 0.7) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~species,scales='fixed', space = "free") +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
  stat_summary(fun.y="mean", geom="point", size=5, 
               aes(shape = physiology,color=physiology)) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size= 20),
          axis.text.x = element_text(angle = 90, hjust = 1,size=20),
          axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 12)) +
    ggtitle("Clade 2")
plot(C2)
C3<-ggplot(tcounts %>% filter(clade=='Clade3'), aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
  ylim(-2.4, 10) +
    geom_point(size=0.7) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~species,scales='fixed', space = "free") +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="errorbar", aes(color=physiology), width=0.2) +
  stat_summary(fun.y="mean", geom="point", size=5, 
               aes(shape = physiology,color=physiology)) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
          axis.text.y = element_text(size= 20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 12)) +
    ggtitle("Clade 3")
plot(C3)
grid.arrange(C1,C2,C3,ncol=3, widths = c(7,3,4))

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
          theme(axis.text.x=element_text(size=16,angle=90, hjust=1),axis.text.y=element_text(size=16),axis.title.y = element_text(size=20),legend.text=element_text(size=20),strip.text.x = element_text(size = 10)))
}
#dev.off()

ext = element_text(size=20),
axis.text.x = element_text(angle=90, hjust=1)

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

goi<-c("ENSFHEP00000017866","ENSFHEP00000033132")
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
# test gene plots
# use this to fill in gene ID and see what it looks like
# -------------------
# OCLN
goi<-c("ENSFHEP00000018175","ENSFHEP00000018194","ENSFHEP00000032878")
# slc
# clcn
goi<-c("ENSFHEP00000019608")
# SLC
#goi<-c("ENSFHEP00000019613","ENSFHEP00000019835","ENSFHEP00000019926","ENSFHEP00000020009","ENSFHEP00000020346","ENSFHEP00000020382","ENSFHEP00000020700","ENSFHEP00000020976","ENSFHEP00000020992","ENSFHEP00000020992","ENSFHEP00000021656","ENSFHEP00000022018","ENSFHEP00000022096","ENSFHEP00000022161","ENSFHEP00000022287","ENSFHEP00000022366","ENSFHEP00000022435","ENSFHEP00000022454","ENSFHEP00000022475","ENSFHEP00000022598","ENSFHEP00000022613","ENSFHEP00000023404","ENSFHEP00000023470","ENSFHEP00000023557","ENSFHEP00000023589","ENSFHEP00000023609")

goi<-c("ENSFHEP00000006227","ENSFHEP00000026135","ENSFHEP00000006110","ENSFHEP00000009214","ENSFHEP00000004076",
       "ENSFHEP00000023773","ENSFHEP00000023784","ENSFHEP00000033997","ENSFHEP00000035103","ENSFHEP00000018763",
       "ENSFHEP00000007353","ENSFHEP00000007376","ENSFHEP00000025947","ENSFHEP00000010833","ENSFHEP00000022487",
       "ENSFHEP00000003286","ENSFHEP00000002741","ENSFHEP00000009396","ENSFHEP00000009409","ENSFHEP00000006579",
       "ENSFHEP00000014084","ENSFHEP00000000720","ENSFHEP00000009893","ENSFHEP00000006308","ENSFHEP00000000898",
       "ENSFHEP00000031738","ENSFHEP00000022435","ENSFHEP00000001295","ENSFHEP00000004102","ENSFHEP00000004089",
       "ENSFHEP00000005954","ENSFHEP00000014786","ENSFHEP00000034578","ENSFHEP00000009999","ENSFHEP00000023609",
       "ENSFHEP00000022287","ENSFHEP00000000795","ENSFHEP00000016688","ENSFHEP00000008627","ENSFHEP00000034356",
       "ENSFHEP00000014782","ENSFHEP00000003491","ENSFHEP00000000933","ENSFHEP00000000615","ENSFHEP00000024212",
       "ENSFHEP00000000292","ENSFHEP00000006166","ENSFHEP00000034768","ENSFHEP00000030251","ENSFHEP00000016967",
       "ENSFHEP00000019926","ENSFHEP00000005103","ENSFHEP00000026744","ENSFHEP00000031637","ENSFHEP00000017152",
       "ENSFHEP00000032431","ENSFHEP00000013059","ENSFHEP00000004785","ENSFHEP00000029010","ENSFHEP00000002249",
       "ENSFHEP00000001350","ENSFHEP00000003852","ENSFHEP00000024546","ENSFHEP00000022475","ENSFHEP00000031723",
       "ENSFHEP00000018441","ENSFHEP00000003611","ENSFHEP00000000608","ENSFHEP00000003596","ENSFHEP00000023850",
       "ENSFHEP00000010805","ENSFHEP00000005530","ENSFHEP00000003883","ENSFHEP00000022161","ENSFHEP00000022325",
       "ENSFHEP00000024569","ENSFHEP00000017403","ENSFHEP00000008715","ENSFHEP00000004556","ENSFHEP00000003724",
       "ENSFHEP00000022598","ENSFHEP00000018824","ENSFHEP00000018834","ENSFHEP00000032132","ENSFHEP00000017082",
       "ENSFHEP00000005485","ENSFHEP00000005503","ENSFHEP00000013010","ENSFHEP00000029651","ENSFHEP00000020009",
       "ENSFHEP00000031934","ENSFHEP00000034198","ENSFHEP00000014052","ENSFHEP00000016137","ENSFHEP00000031733",
       "ENSFHEP00000020346","ENSFHEP00000028832","ENSFHEP00000004993","ENSFHEP00000033204","ENSFHEP00000029959",
       "ENSFHEP00000020382","ENSFHEP00000026760","ENSFHEP00000031658","ENSFHEP00000019932","ENSFHEP00000014036",
       "ENSFHEP00000028223","ENSFHEP00000000067","ENSFHEP00000032217","ENSFHEP00000012529","ENSFHEP00000004837",
       "ENSFHEP00000010074","ENSFHEP00000027737","ENSFHEP00000027791","ENSFHEP00000013771","ENSFHEP00000030175",
       "ENSFHEP00000030191","ENSFHEP00000006780","ENSFHEP00000010150","ENSFHEP00000003964","ENSFHEP00000018233",
       "ENSFHEP00000018241","ENSFHEP00000004473","ENSFHEP00000005872","ENSFHEP00000002602","ENSFHEP00000018374",
       "ENSFHEP00000018387","ENSFHEP00000018399","ENSFHEP00000000821","ENSFHEP00000014149","ENSFHEP00000001072",
       "ENSFHEP00000024125","ENSFHEP00000023557","ENSFHEP00000033959","ENSFHEP00000020104","ENSFHEP00000022750",
       "ENSFHEP00000022454","ENSFHEP00000004169","ENSFHEP00000022366","ENSFHEP00000007397","ENSFHEP00000002304",
       "ENSFHEP00000022613","ENSFHEP00000004992","ENSFHEP00000030432","ENSFHEP00000020992",'ENSFHEP00000033559',
       "ENSFHEP00000015049","ENSFHEP00000031020","ENSFHEP00000001078","ENSFHEP00000001090","ENSFHEP00000015367",
       "ENSFHEP00000031245","ENSFHEP00000012623","ENSFHEP00000022096","ENSFHEP00000000143","ENSFHEP00000002610",
       "ENSFHEP00000005895","ENSFHEP00000010032","ENSFHEP00000029267","ENSFHEP00000000330","ENSFHEP00000016685",
       "ENSFHEP00000016699","ENSFHEP00000016718","ENSFHEP00000028108","ENSFHEP00000015356","ENSFHEP00000028837",
       "ENSFHEP00000026805","ENSFHEP00000000800","ENSFHEP00000020700","ENSFHEP00000006237","ENSFHEP00000008921",
       "ENSFHEP00000008929","ENSFHEP00000020976","ENSFHEP00000029859","ENSFHEP00000002444","ENSFHEP00000018950",
       "ENSFHEP00000013896","ENSFHEP00000012482","ENSFHEP00000029285","ENSFHEP00000005913","ENSFHEP00000024967",
       "ENSFHEP00000029789","ENSFHEP00000014508","ENSFHEP00000013200","ENSFHEP00000019613","ENSFHEP00000034022",
       "ENSFHEP00000006989","ENSFHEP00000004489","ENSFHEP00000001354","ENSFHEP00000019359","ENSFHEP00000022541",
       "ENSFHEP00000025379","ENSFHEP00000015425","ENSFHEP00000028455","ENSFHEP00000032648","ENSFHEP00000011768",
       "ENSFHEP00000019835","ENSFHEP00000026156","ENSFHEP00000010040","ENSFHEP00000024007","ENSFHEP00000017788",
       "ENSFHEP00000027123","ENSFHEP00000023752","ENSFHEP00000027097","ENSFHEP00000014068","ENSFHEP00000030094",
       "ENSFHEP00000013674","ENSFHEP00000013373","ENSFHEP00000020051","ENSFHEP00000016313","ENSFHEP00000013140",
       "ENSFHEP00000005147","ENSFHEP0000000516","ENSFHEP00000030091","ENSFHEP00000030128","ENSFHEP00000026043",
       "ENSFHEP00000021656","ENSFHEP00000013385","ENSFHEP00000026601","ENSFHEP00000011422","ENSFHEP00000022882",
       "ENSFHEP00000005888","ENSFHEP00000005702","ENSFHEP00000024852","ENSFHEP00000029030","ENSFHEP00000030417",
       "ENSFHEP00000030425","ENSFHEP00000005899","ENSFHEP00000005968","ENSFHEP00000006035","ENSFHEP0000003475",
       "ENSFHEP00000022912","ENSFHEP00000000129","ENSFHEP00000025237","ENSFHEP00000013790","ENSFHEP00000030070")
#goi <- c("ENSFHEP00000015765","ENSFHEP00000004305","ENSFHEP00000033970","ENSFHEP00000016264","ENSFHEP00000016277","ENSFHEP00000017105","ENSFHEP00000032374","ENSFHEP00000030512","ENSFHEP00000007804","ENSFHEP00000026234","ENSFHEP00000019144")
# three-way
threeway <- rownames(sig_threeway)
goi <- threeway
goi <- c("ENSFHEP00000006340")
goi <- c("ENSFHEP00000004188")
goi <- c("ENSFHEP00000023858")
goi <- c("ENSFHEP00000005982")
goi <- c("ENSFHEP00000018189")
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

  