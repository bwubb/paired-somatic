#!/usr/bin/env Rscript
#Merge per-chr CODEX2 filtered segment tables.
#
#Rscript codex2_merge.R --outdir data/work/codex2 --output data/work/codex2/codex2.segments.filtered.txt

parse_args<-function(argv){
    out<-list()
    i<-1L
    while(i<=length(argv)){
        a<-argv[i]
        if(!startsWith(a,"--")) stop("expected --flag, got: ",a)
        key<-substring(a,3L)
        if(i==length(argv)||startsWith(argv[i+1L],"--")){
            out[[key]]<-TRUE
            i<-i+1L
        }else{
            out[[key]]<-argv[i+1L]
            i<-i+2L
        }
    }
    out
}

args<-parse_args(commandArgs(trailingOnly=TRUE))
if(is.null(args$outdir)||is.null(args$output)) stop("Required: --outdir --output")
outdir<-args$outdir
outfile<-args$output

files<-sort(list.files(outdir,pattern="^chr.+\\.codex2\\.segments\\.filtered\\.txt$",full.names=TRUE))
if(!length(files)) stop("no filtered segment files in ",outdir)

tabs<-lapply(files,function(f){
    x<-tryCatch(read.delim(f,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
    if(is.null(x)||!nrow(x)) return(NULL)
    x
})
tabs<-Filter(Negate(is.null),tabs)
if(!length(tabs)){
    write.table(data.frame(),file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
    message("No CNV rows after filter; wrote empty ",outfile)
    quit(save="no",status=0)
}
out<-do.call(rbind,tabs)
write.table(out,file=outfile,sep="\t",quote=FALSE,row.names=FALSE)
message("Merged ",length(tabs)," chr tables -> ",outfile," (",nrow(out)," rows)")
