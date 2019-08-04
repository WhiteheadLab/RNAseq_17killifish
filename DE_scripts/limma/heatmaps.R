library(pheatmap)
library(dplyr)

######## IMPORT DATA ###########

setwd("~/Documents/Whitehead_Office-Desktop/Research/NSF-osmotic_tolerance/15.species/Lisa.Johnson.data")
df <- data.frame(read.csv("normalized_counts_V2.csv", row.name=1))

######## SELECT DATA USING PVALUES, AND SELECT COLUMNS TO PLOT ###########

mydata_pvalueselect <- filter(df, Salin.Physiol.adjP<=0.05)
#mydata_pvalueselect <- filter(df, Salin.Physiol.adjP<=0.01 & X3way.adjP>0.1)
#mydata_pvalueselect <- filter(df, Salinity.adjP<=0.01 & Salin.Physiol.adjP>0.1 & X3way.adjP>0.1)
#mydata_pvalueselect <- filter(df, Salinity.adjP<=0.01 & Salin.Physiol.adjP>0.1 & X3way.adjP>0.1 & Average_pos.neg=="pos")

mydata_log2expr_relative_to_15ppt_per_species <- data.frame(mydata_pvalueselect[,58:85])
mydata1 <- as.matrix(mydata_log2expr_relative_to_15ppt_per_species)

######## TREATMENT LABELING ###########

sample <- c("Fsim.M.0","Fsim.M.15","Fgds.M.0","Fgds.M.15","FheL.M.0","FheL.M.15","FheP.M.0","FheP.M.15","Fdia.M.0","Fdia.M.15","Fcat.F.0","Fcat.F.15","Frat.F.0","Frat.F.15","Fpns.M.0","Fpns.M.15","Lpva.M.0","Lpva.M.15","Lgdi.F.0","Lgdi.F.15","Fxen.M.0","Fxen.M.15","Fchr.M.0","Fchr.M.15","Fnts.F.0","Fnts.F.15","Foli.F.0","Foli.F.15")
clade <- c("clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade1","clade2","clade2","clade2","clade2","clade2","clade2","clade3","clade3","clade3","clade3","clade3","clade3","clade3","clade3")
physiology <- c("M","M","M","M","M","M","M","M","M","M","F","F","F","F","M","M","M","M","F","F","M","M","M","M","F","F","F","F")
salinity <- c("0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt","0ppt","15ppt")
sample_label_df <- data.frame(salinity,physiology,clade)
#add sample_label as the row names for the above data.frame
.rowNamesDF(sample_label_df, make.names=FALSE) <- sample

annotation_colors = list(salinity = c("0ppt"="red", "15ppt"="blue"),physiology = c("M"="blue", "F"="red"),clade = c("clade1"="green", "clade2"="orange", "clade3"="purple"))

######## SET HEATMAP COLOR SCALE ###########

my.breaks <- c(seq(-2, 0, by=0.1), seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("#00FFFF", "black"))(length(my.breaks)/2), colorRampPalette(colors = c("black", "yellow"))(length(my.breaks)/2))

######## GENERATE HEATMAP ###########

pheatmap(mydata1, cluster_rows = TRUE,
         clustering_distance_rows = "correlation",
         cluster_cols = FALSE,
         annotation_col = sample_label_df,
         show_rownames=F,
         cutree_rows = 3,
         color = my.colors,
         annotation_colors = annotation_colors,
         breaks = my.breaks,
         gaps_col = c(14,20),
         border_color = NA,
)