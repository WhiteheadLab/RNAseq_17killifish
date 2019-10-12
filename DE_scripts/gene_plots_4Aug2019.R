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


# save for when you want to summarize by species_condition
#species_condition<-as.vector(paste(species,condition,sep="."))
#colnames(norm_counts) <- species_condition
#dim(norm_counts) 
#mean_norm_counts<-sapply(unique(colnames(norm_counts)), function(i)
#  rowMeans(norm_counts[,colnames(norm_counts) == i]))

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


# ----------------
# single genes plots
# ------------------

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


# ENSFHEP00000033747
# per3
goi <- c("ENSFHEP00000019155")
# use this one:
goi <- c("ENSFHEP00000033747")
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

# ENSFHEP00000006725
# aqp3
goi<-c("ENSFHEP00000006725")
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










# circadian genes
circ_id <- unique(ann[ann$ensembl_gene_id %in% circ_genes,]$ensembl_peptide_id)
circ_id <- circ_id[!circ_id == "ENSFHEP00000019155"]
goi <- circ_id
# per1B
goi <- c("ENSFHEP00000023694")
# nr1d1
goi <- c("ENSFHEP00000032821")
# bhlhe40
goi <- c("ENSFHEP00000032833")
# per3
goi <- c("ENSFHEP00000033747")

# physiology genes
#ext
goi <- c("ENSFHEP00000018189")
#slc25a12
goi <- c("ENSFHEP00000017152")
# vti1a
goi <- c("ENSFHEP00000030532")
# exoc7
goi <- c("ENSFHEP00000009514")
# slc35a3a
goi <- c("ENSFHEP00000001090")

goi <- c("ENSFHEP00000008352")

# taok1b
goi <- c("ENSFHEP00000028829")
# atp1a1a
goi <- c("ENSFHEP00000015386")


tmp <- norm_counts[rownames(norm_counts) %in% goi,]
tmp_ann <- merge(tmp,ann,by.x = "row.names", by.y  = "ensembl_peptide_id")
dim(tmp_ann)
tmp_ann_tmp <- tmp_ann[!duplicated(tmp_ann$Row.names),]
rownames(tmp_ann_tmp) <- tmp_ann_tmp$external_gene_name



tmp <- tmp_ann[,c(2:82)]
rownames(tmp) <- tmp_ann_tmp$external_gene_name
tmp <- data.frame(tmp,stringsAsFactors = FALSE)
tmp <- data.matrix(tmp)
head(tmp,quote = FALSE)
ltmp <- log2(tmp+0.5)
tcounts <- t(ltmp) %>%
  merge(ExpDesign, ., by="row.names") %>% 
  gather(gene, expression, (ncol(.)-length(rownames(tmp))+1):ncol(.))
tcounts %>% dplyr::select(Row.names, clade, physiology, condition, gene, expression) %>% head %>% knitr::kable() %>% kable_styling()


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




# find NKCC
# NKATP a1
# NKATP a2

