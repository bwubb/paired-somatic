#!/usr/bin/env Rscript
#Germline CODEX2 per-chromosome: normalize_codex2_ns + integer CBS.
#
#Rscript codex2_run_chr.R \
#  --outdir data/work/codex2 \
#  --chr 13 \
#  --controls data/work/codex2/controls.list \
#  [--Kmax 10]

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
if(is.null(args$outdir)||is.null(args$chr)||is.null(args$controls)){
    stop("Required: --outdir --chr --controls")
}
outdir<-args$outdir
chr<-args$chr
controls_fp<-args$controls
Kmax<-if(!is.null(args$Kmax)) as.integer(args$Kmax) else 10L
K<-seq_len(Kmax)

suppressPackageStartupMessages(library(CODEX2))

need<-c("coverageQC.csv","ref_qc.rds","library_size_factor.csv")
for(f in need){
    fp<-file.path(outdir,f)
    if(!file.exists(fp)) stop("missing prep output: ",fp)
}
if(!file.exists(controls_fp)) stop("controls list not found: ",controls_fp)

Y_qc<-as.matrix(read.csv(file.path(outdir,"coverageQC.csv"),row.names=1,check.names=FALSE))
ref_qc<-readRDS(file.path(outdir,"ref_qc.rds"))
Ndf<-read.csv(file.path(outdir,"library_size_factor.csv"),check.names=FALSE)
if("sample"%in%names(Ndf)&&"N"%in%names(Ndf)){
    N<-setNames(as.numeric(Ndf$N),as.character(Ndf$sample))
    N<-N[colnames(Y_qc)]
}else{
    N<-as.numeric(Ndf[,ncol(Ndf)])
    names(N)<-colnames(Y_qc)
}
if(any(is.na(N))) stop("library size N does not align to coverageQC columns")

controls<-scan(controls_fp,what=character(),quiet=TRUE)
controls<-controls[nzchar(controls)&!grepl("^#",controls)]
norm_index<-which(colnames(Y_qc)%in%controls)
if(!length(norm_index)) stop("no control samples found in coverageQC colnames")
missing_ctrl<-setdiff(controls,colnames(Y_qc))
if(length(missing_ctrl)){
    warning("controls not in QC matrix (dropped earlier?): ",paste(missing_ctrl,collapse=","))
}
message("chr=",chr," exons=",length(ref_qc)," samples=",ncol(Y_qc)," controls_in_norm=",length(norm_index))

chr.index<-which(as.character(seqnames(ref_qc))==as.character(chr))
if(!length(chr.index)) stop("no QC exons for chr ",chr)

sampname_qc<-colnames(Y_qc)
gc_qc<-ref_qc$gc
if(is.null(gc_qc)) stop("ref_qc lacks gc metadata; re-run prep")

message("normalize_codex2_ns K=1:",Kmax," ...")
normObj<-normalize_codex2_ns(
    Y_qc=Y_qc[chr.index,,drop=FALSE],
    gc_qc=gc_qc[chr.index],
    K=K,
    norm_index=norm_index,
    N=N
)
Yhat.ns<-normObj$Yhat
BIC.ns<-normObj$BIC
AIC.ns<-normObj$AIC
RSS.ns<-normObj$RSS

choiceofK(AIC.ns,BIC.ns,RSS.ns,K=K,filename=file.path(outdir,paste0("codex2_chr",chr,"_ns_choiceofK.pdf")))

optK<-which.max(BIC.ns)
message("optK=",optK," (max BIC)")

finalcall.CBS<-segmentCBS(
    Y_qc[chr.index,,drop=FALSE],
    Yhat.ns,
    optK=optK,
    K=K,
    sampname_qc=sampname_qc,
    ref_qc=ranges(ref_qc)[chr.index],
    chr=chr,
    lmax=400,
    mode="integer"
)
out_raw<-file.path(outdir,paste0("chr",chr,".codex2.segments.txt"))
write.table(finalcall.CBS,file=out_raw,sep="\t",quote=FALSE,row.names=FALSE)

if(nrow(finalcall.CBS)){
    filter1<-finalcall.CBS$length_kb<=200
    filter2<-finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
    finalcall.CBS.filter<-finalcall.CBS[filter1&filter2,,drop=FALSE]
    if(nrow(finalcall.CBS.filter)){
        filter3<-finalcall.CBS.filter$lratio>40
        filter4<-(finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
        finalcall.CBS.filter<-finalcall.CBS.filter[filter3|filter4,,drop=FALSE]
    }
}else{
    finalcall.CBS.filter<-finalcall.CBS
}
out_filt<-file.path(outdir,paste0("chr",chr,".codex2.segments.filtered.txt"))
write.table(finalcall.CBS.filter,file=out_filt,sep="\t",quote=FALSE,row.names=FALSE)
message("Wrote ",out_raw," and ",out_filt," (n_raw=",nrow(finalcall.CBS)," n_filt=",nrow(finalcall.CBS.filter),")")
