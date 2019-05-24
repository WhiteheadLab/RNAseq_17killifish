library(gridExtra)
library(data.table)
library(dplyr)
library(grid)
library(ggplot2)
library(lattice)
library(biomaRt)


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

# contrasts

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
# main and interaction effects, F test
main_phys <- read.csv("main_phys.csv", stringsAsFactors = FALSE, header = TRUE)
main_clade <- read.csv("main_clade.csv", stringsAsFactors = FALSE, header = TRUE)
main_condition <- read.csv("main_condition.csv", stringsAsFactors = FALSE, header = TRUE)
int_phys_condition <- read.csv("interaction_physiology_condition.csv", stringsAsFactors = FALSE, header = TRUE)
int_phys_clade <- read.csv("interaction_physiology_clade.csv", stringsAsFactors = FALSE, header = TRUE)
int_clade_condition <- read.csv("interaction_clade_condition.csv", stringsAsFactors = FALSE, header = TRUE)
int_three <- read.csv("interaction_threeway.csv", stringsAsFactors = FALSE, header = TRUE)

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

# significant genes
# main effects
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

d<-as.matrix(mean_norm_counts_ordered_clade_sig)
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

rld <- log2(mean_norm_counts_ordered_clade_sig+1)
geneDists <- dist(mean_norm_counts_ordered_clade_sig)
df <- data.frame(ph,cl, condition,stringsAsFactors=FALSE)
rownames(df) <- colnames(rld)
pheatmap(rld, show_rownames=FALSE,
         clustering_distance_rows = geneDists, 
         cluster_cols= FALSE,
         annotation_col=df,
         scale = "none")
head(rld)
# for each FW col, subtract values by itself and subtract 15ppt col by original FW
# in for loop

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

all_goi<-c("ENSFHEP00000007220.1","ENSFHEP00000025841","ENSFHEP00000019510",
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
  
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='fixed') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='fixed') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='fixed',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}


# ----------------------------
# main effects
length(sig_main_phys)
tmp <- j[rownames(j) %in% sig_main_phys,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
#  geom_boxplot() + 
#  facet_wrap(~gene,scales="free_x") +
#  labs(x="Clade",
#       y="Expression of sig for three-way (log2 cpm normalized counts)",
#       fill = "Native Physiology") +
#  theme_bw()


pdf("~/Documents/UCDavis/Whitehead/physiology_maineffect_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()

# condition
length(sig_main_condition)
tmp <- j[rownames(j) %in% sig_main_condition,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene,scales="free_x") +
  labs(x="Clade",
       y="Expression of sig for three-way (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw()


pdf("~/Documents/UCDavis/Whitehead/condition_maineffect_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()


# clade
length(sig_main_clade)
# [1] 2614
tmp <- j[rownames(j) %in% sig_main_clade,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
#  geom_boxplot() + 
#  facet_wrap(~gene,scales="free_x") +
#  labs(x="Clade",
#       y="Expression of sig for three-way (log2 cpm normalized counts)",
#       fill = "Native Physiology") +
#  theme_bw()

#ggplot(tcounts,aes(x=clade, y=expression,fill=physiology)) + 
#  geom_boxplot() + 
#  facet_wrap(~gene,scales="free_x") +
#  labs(x="Clade",
#       y="Expression of sig for three-way (log2 cpm normalized counts)",
#       fill = "Native Physiology") +
#  theme_bw()

pdf("~/Documents/UCDavis/Whitehead/clade_maineffect_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()

# ---------------------------
# two-way interactions
length(sig_int_phys_condition)
tmp <- j[rownames(j) %in% sig_int_phys_condition,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene,scales="free_x") +
  labs(x="Clade",
       y="Expression of sig for three-way (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw()

pdf("~/Documents/UCDavis/Whitehead/physiologycondition_twowayinteraction_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()


length(sig_int_phys_clade)
tmp <- j[rownames(j) %in% sig_int_phys_clade,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

#ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
#  geom_boxplot() + 
#  facet_wrap(~gene,scales="free_x") +
#  labs(x="Clade",
#       y="Expression of sig for three-way (log2 cpm normalized counts)",
#       fill = "Native Physiology") +
#  theme_bw()


pdf("~/Documents/UCDavis/Whitehead/physiologyclade_twowayinteraction_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()




length(sig_int_clade_condition)
tmp <- j[rownames(j) %in% sig_int_clade_condition,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,alpha=clade,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene,scales="free_x") +
  labs(y="Expression (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw()

pdf("~/Documents/UCDavis/Whitehead/cladecondition_twowayinteraction_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()


# --------------------------
# three-way interaction
length(sig_int_three)
sig_int_three

tmp <- j[rownames(j) %in% sig_int_three,]
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()

ggplot(tcounts,aes(x=clade:condition, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene,scales="free_x") +
  labs(x="Clade",
       y="Expression (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw()

ggplot(tcounts,aes(x=clade, y=expression,fill=physiology)) + 
  geom_boxplot() + 
  facet_wrap(~gene,scales="free_x") +
  labs(x="Clade",
       y="Expression of sig for three-way (log2 cpm normalized counts)",
       fill = "Native Physiology") +
  theme_bw()

pdf("~/Documents/UCDavis/Whitehead/physiologycladecondition_threewayinteraction_23May2019.pdf",paper="USr",width=13.5, height=8)
for (i in rownames(ltmp)){
  new <- data.frame(ltmp[i, ])
  colnames(new) <- c(i)
  tcounts <- new %>% 
    merge(ExpDesign, ., by="row.names") %>% 
    gather(gene, expression, (ncol(.)-length(colnames(new))+1):ncol(.))
  
  C1<-ggplot(tcounts %>%
               filter(clade=='Clade1'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) +
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(y="Expression (log2 normalized counts)")+
    ggtitle("Clade 1")
  
  #plot(C1)
  C2<-ggplot(tcounts %>%
               filter(clade=='Clade2'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line",aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free') +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology),width=0.2) +
    theme_bw() +
    theme(legend.position="bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          strip.text.y = element_blank()) +
    labs(x="salinity treatment")+
    ggtitle("Clade 2")
  #plot(C2)
  C3<-ggplot(tcounts %>%
               filter(clade=='Clade3'),
             aes(factor(condition,levels = c("0.2_ppt","15_ppt")), expression)) + 
    geom_point(aes(color=physiology)) +
    stat_summary(fun.y="mean", geom="line", aes(group=physiology,color=physiology)) +
    facet_grid(~gene,scales='free',labeller=) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", aes(color=physiology), width=0.2) +
    theme_bw() +
    theme(legend.position="none",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    ggtitle("Clade 3")
  #plot(C3)
  grid.arrange(C1,C2,C3,ncol=3)
}
dev.off()
