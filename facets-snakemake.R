###facets-snakemake.R

library("facets")

#init
datafile=snakemake@input[[1]]
cval=snakemake@params[["cval"]]
outpath=dirname(normalizePath(datafile))

#run
rcmat=readSnpMatrix(datafile)
xx=preProcSample(rcmat)
oo=procSample(xx,cval=150)
fit=emcncf(oo)

#save pdf
pdf(paste(outpath,"copynumber_profile.pdf",sep="/"))
plotSample(x=oo,emfit=fit)
dev.off()
pdf(paste(outpath,"fit_diagnostic.pdf",sep="/"))
logRlogORspider(oo$out,oo$dipLogR)
dev.off()

#write output
write.csv(data.frame(purity=fit$purity,ploidy=fit$ploidy),file=paste(outpath,"purity_ploidy.csv",sep="/"),row.names=FALSE,quote=FALSE)
write.csv(fit$cncf,file=paste(outpath,"segmentation_cncf.csv",sep="/"),row.names=FALSE,quote=FALSE)