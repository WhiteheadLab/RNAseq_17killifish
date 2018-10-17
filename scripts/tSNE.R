# remove F.het.MDPP and F.het.MDPL and F.grandis from PCA
# what are % PCA
library(ggplot2)

design
sp<-as.character(unlist(design[1,]))
sp<-sp[-c(1,2)]
ph<-as.character(unlist(design[2,]))
ph<-ph[-c(1,2)]
cl<-as.character(unlist(design[3,]))
cl<-cl[-c(1,2)]
de<-as.character(unlist(design[4,]))
de<-de[-c(1,2)]
# clade
names<-colnames(log_x)
tsne<-Rtsne(t(log_x),dims=2,perplexity=10,verbose=T,max_iter=1000)
tplot<-cbind(v1=tsne$Y[,1],v2=tsne$Y[,2],clade)
a<-as.data.frame(tplot)
ggplot(a,
       aes(x=v1,y=v2,color=cl,label=names))+
  geom_point(cex=3) +
  geom_text(aes(label=names),hjust=0,vjust=2)+
  theme_classic() +
  labs(x="tsne1",y="tsne2")+
  theme(axis.line=element_line(size=1.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text=element_text(size=15))
# genes
set.seed(5)
tsne<-Rtsne(log_x,dims=2,perplexity=50,verbose=T,max_iter=1000,check_duplicates = FALSE)
tplot<-cbind(v1=tsne$Y[,1],v2=tsne$Y[,2])
a<-as.data.frame(tplot)
ggplot(a,
       aes(x=v1,y=v2))+
  geom_point(cex=1) +
  theme_classic() +
  labs(x="tsne1",y="tsne2")+
  theme(axis.line=element_line(size=1.5),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text=element_text(size=15))
