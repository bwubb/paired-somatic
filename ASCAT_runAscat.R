library(ASCAT)
#Note: They switched back to Tumor...

ascat.bc = ascat.loadData(
    Tumor_LogR_file=snakemake@input[["tumorLogR"]],
    Tumor_BAF_file=snakemake@input[["tumorBAF"]],
    Germline_LogR_file=snakemake@input[["germlineLogR"]],
    Germline_BAF_file=snakemake@input[["germlineBAF"]],
    gender='XX',
    genomeVersion="hg38",
    isTargetedSeq=T)

ascat.plotRawData(ascat.bc,img.dir=snakemake@params[["outdir"]],img.prefix="Before_correction_")
ascat.bc=ascat.correctLogR(ascat.bc,GCcontentfile = "GC_G1000_hg38.txt")
ascat.plotRawData(ascat.bc,img.dir=snakemake@params[["outdir"]],img.prefix="After_correction_")
ascat.bc=ascat.aspcf(ascat.bc,out.dir=snakemake@params[["outdir"]],penalty=25)
ascat.plotSegmentedData(ascat.bc,img.dir=snakemake@params[["outdir"]])
ascat.output=ascat.runAscat(ascat.bc,gamma=1,write_segments=T,img.dir=snakemake@params[["outdir"]])
QC=ascat.metrics(ascat.bc,ascat.output)
write.table(QC,file=paste(c(snakemake@params[["outdir"]],"QC.txt"),sep="/"),sep=" ")
save(ascat.bc,ascat.output,QC,file=snakemake@output[["rdata"]])
