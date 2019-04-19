rld <- vst(dds, blind = FALSE,fitType='local')
sampleDists <- dist(t(assay(rld)))
df <- as.data.frame(colData(dds)[,c("physiology","condition","clade")])
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, annotation = df, show_rownames=F)

select100 <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]

pheatmap(assay(rld)[select100,], show_rownames=T,clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, annotation_col=df)

select500 <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:500]

sampleDists <- dist(t(assay(rld)[select500,]))
sampleDistMatrix <- as.matrix( sampleDists )

pheatmap(sampleDistMatrix, show_rownames=T,clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, annotation_col=df)



cowplot::plot_grid( plotPCA(rld, intgroup="physiology"),
                    plotPCA(rld, intgroup="condition"),
                    plotPCA(rld, intgroup=c( "clade")),
                    align="c", ncol=2)
