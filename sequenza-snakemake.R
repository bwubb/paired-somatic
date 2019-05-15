############################
#Brad Wubbenhorst
#bwubb@pennmedicine.upenn.edu
#Dec. 2018

library(sequenza)
#data.file <- system.file("data", "example.seqz.txt.gz", package = "sequenza")
tumor<-sequenza.extract(snakemake@input[[1]])
CP<-sequenza.fit(tumor,mc.cores=snakemake@threads)
sequenza.results(sequenza.extract=tumor,cp.table=CP,sample.id=snakemake@wildcards[["tumor"]],out.dir=snakemake@params[["outdir"]])