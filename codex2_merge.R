#!/usr/bin/env Rscript
#Merge per-chr CODEX2 filtered segment tables.
#Usage: Rscript codex2_merge.R <outdir> <outfile>

args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("Usage: Rscript codex2_merge.R <outdir> <outfile>")
outdir<-args[1]
outfile<-args[2]

files<-sort(list.files(outdir,pattern="^chr.+\\.codex2\\.segments\\.filtered\\.txt$",full.names=TRUE))
if(!length(files)) stop("no filtered segment files in ",outdir)

tabs<-lapply(files,function(f){
    x<-tryCatch(read.delim(f,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
    if(is.null(x)||!nrow(x)) return(NULL)
    x
})
tabs<-Filter(Negate(is.null),tabs)
if(!length(tabs)){
    #Still write empty header-ish file for Snakemake.
    write.table(data.frame(),file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
    message("No CNV rows after filter; wrote empty ",outfile)
    quit(save="no",status=0)
}
out<-do.call(rbind,tabs)
write.table(out,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
message("Merged ",length(tabs)," chr tables -> ",outfile," (",nrow(out)," rows)")
