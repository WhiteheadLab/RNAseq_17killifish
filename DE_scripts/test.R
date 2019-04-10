# ---------------------
# dream tutorial:
# https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html
# this is for repeated measures designs (where same patient tested twice for same variable, like control measurement first then treatment)
# ---------------------

library('variancePartition')
library('edgeR')
library('doParallel')
data(varPartDEdata)
isexpr = rowSums(cpm(countMatrix)>0.1) >= 3
genes = DGEList( countMatrix )
genes = calcNormFactors( genes )
design = model.matrix( ~ Disease, metadata)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata$Individual)
vobj = voom( genes, design, plot=FALSE, block=metadata$Individual, correlation=dupcor$consensus)
design = model.matrix( ~ Disease, metadata)
dupcor <- duplicateCorrelation(vobj, design, block=metadata$Individual)
fitDupCor <- lmFit(vobj, design, block=metadata$Individual, correlation=dupcor$consensus)
fitDupCor <- eBayes( fitDupCor )
cl <- makeCluster(4)
registerDoParallel(cl)
form <- ~ Disease + (1|Individual) 
L = getContrast( vobj, form, metadata, "Disease1")
# took 274 s on 4/10/2019
fitmm = dream( vobj, form, metadata, L)
fitmm = eBayes( fitmm )
topTable( fitmm )
form <- ~ 0 + Disease + (1|Individual) 
L = getContrast( vobj, form, metadata, c("Disease1", "Disease0"))
L
form <- ~ 0 + Disease + (1|Individual) 

# define and then cbind contrasts
L1 = getContrast( vobj, form, metadata, "Disease0")
L2 = getContrast( vobj, form, metadata, "Disease1")
L = cbind(L1, L2)

# fit both contrasts
fit = dream( vobj[1:10,], form, metadata, L)
fiteb = eBayes( fit )

# extract results from first contrast
topTable( fiteb, coef="L1" )
form = ~ (1|Individual) + (1|Disease)
vp = fitExtractVarPartModel( vobj, form, metadata)
plotVarPart( sortCols(vp))


# Compare p-values and make plot
p1 = topTable(fitDupCor, coef="Disease1", number=Inf, sort.by="none")$P.Value
p2 = topTable(fitmm, number=Inf, sort.by="none")$P.Value

plotCompareP( p1, p2, vp$Individual, dupcor$consensus)
