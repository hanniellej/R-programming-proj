library("affy")
library("limma")
pd = read.AnnotatedDataFrame("samples.txt",header=TRUE,sep=",",row.names=1)
mydata = ReadAffy(filenames=pd$filename,phenoData=pd,verbose=TRUE)
sampleNames(mydata)=row.names(pData(pd))
eset = rma(mydata)
factor.vector=factor(c('m','m','m','d','d','d'))
factor.vector=relevel(factor.vector,ref='m')
exp_design=model.matrix(~0 + factor.vector) # formula of linear model that tells R 
colnames(exp_design)=levels(factor.vector)
fit=lmFit(eset,exp_design)
cont.matrix= makeContrasts(
  dvsm = d-m, # difference between mock treated & dicamba treated 
  levels = exp_design)
fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)
N=dim(eset)[1]
results=topTable(fit2,coef=1,number=N,p.value=1)
gene_up=results[results$logFC>=1,]
gene_up_sig=sum(gene_up$adj.P.Value<=0.05)
gene_down=results[results$logFC<=(-1),] # gene that are downregulated
gene_down_sig=sum(gene_down$adj.P.Val<=-0.05)

color="red"

library(RColorBrewer)
barplot(gene_up$AveExpr, beside=TRUE, col="red", main= "Average Differential Expression of Up Regulating Genes", xlab = "Genes", ylab = "Average Expression")
barplot(gene_up$AveExpr,col="green", main= "Average Differential Expression of Down Regulating Genes", xlab = "Genes", ylab = "Average Expression")