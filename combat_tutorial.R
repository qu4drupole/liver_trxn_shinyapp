library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

pheno = pData(bladderEset)
edata = exprs(bladderEset)
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
