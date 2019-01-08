############################
#Brad Wubbenhorst
#bwubb@pennmedicine.upenn.edu
#Dec. 2018

library(sequenza)
tumor<-sequenza.extract(snakemake@input[[1]])
CP<-sequenza.fit(tumor,mc.cores=snakemake@threads)
sequenza.results(sequenza.extract=tumor,cp.table=CP,sample.id=snakemake@wildcards[["tumor"]],out.dir=snakemake@params[["outdir"]])