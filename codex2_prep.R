#!/usr/bin/env Rscript
#Germline CODEX2 prep (GRCh38, contig names without "chr").
#Coverage + QC + library size. Normals only (no tumors).
#
#Usage:
#  Rscript codex2_prep.R <bam.list> <bed> <outdir> [mapq] [mapp]
#    mapq default 20
#    mapp: ones (default) | getmapp

args<-commandArgs(trailingOnly=TRUE)
if(length(args)<3){
    stop("Usage: Rscript codex2_prep.R <bam.list> <bed> <outdir> [mapq] [mapp]")
}
bam_list<-args[1]
bedFile<-args[2]
outdir<-args[3]
mapq<-if(length(args)>=4) as.integer(args[4]) else 20L
mapp_mode<-if(length(args)>=5) args[5] else "ones"

suppressPackageStartupMessages({
    library(CODEX2)
    library(BSgenome.Hsapiens.NCBI.GRCh38)
})

if(!file.exists(bam_list)) stop("bam list not found: ",bam_list)
if(!file.exists(bedFile)) stop("bed not found: ",bedFile)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

bamdir<-scan(bam_list,what=character(),quiet=TRUE)
bamdir<-bamdir[nzchar(bamdir)&!grepl("^#",bamdir)]
if(!length(bamdir)) stop("empty bam list")
missing<-bamdir[!file.exists(bamdir)]
if(length(missing)) stop("missing BAM(s):\n",paste(missing,collapse="\n"))

sampname<-basename(bamdir)
sampname<-sub("\\.ready\\.bam$","",sampname)
sampname<-sub("\\.bam$","",sampname)
if(anyDuplicated(sampname)) stop("duplicate sample names after basename strip")

message("Samples: ",length(sampname))
message("BED: ",bedFile)
message("Outdir: ",outdir)

genome<-BSgenome.Hsapiens.NCBI.GRCh38
bambedObj<-getbambed(bamdir=bamdir,bedFile=bedFile,sampname=sampname,projectname=outdir)
bamdir<-bambedObj$bamdir
sampname<-bambedObj$sampname
ref<-bambedObj$ref
projectname<-bambedObj$projectname

chr_like<-grep("^chr",as.character(seqlevels(ref)),value=TRUE)
if(length(chr_like)){
    stop("BED/ref has chr* contigs (e.g. ",chr_like[1],"). Use GRCh38 BED without chr prefix.")
}

message("Computing GC (NCBI GRCh38)...")
gc<-getgc(ref,genome=genome)
if(identical(mapp_mode,"getmapp")){
    message("Computing mappability via getmapp()...")
    mapp<-tryCatch(getmapp(ref,genome=genome),error=function(e){
        warning("getmapp failed (",conditionMessage(e),"); falling back to mapp=1")
        rep(1,length(ref))
    })
}else{
    message("Setting mappability to 1.")
    mapp<-rep(1,length(ref))
}
values(ref)<-cbind(values(ref),DataFrame(gc,mapp))
saveRDS(ref,file=file.path(projectname,"ref.rds"))

message("Counting coverage (mapq>=",mapq,")...")
coverageObj<-getcoverage(bambedObj,mapqthres=mapq)
Y<-coverageObj$Y
write.csv(Y,file=file.path(projectname,"coverage.csv"),quote=FALSE)

message("QC...")
qcObj<-qc(Y,sampname,ref,cov_thresh=c(20,4000),length_thresh=c(20,2000),mapp_thresh=0.9,gc_thresh=c(20,80))
Y_qc<-qcObj$Y_qc
sampname_qc<-qcObj$sampname_qc
ref_qc<-qcObj$ref_qc
qcmat<-qcObj$qcmat
gc_qc<-ref_qc$gc

write.csv(Y_qc,file=file.path(projectname,"coverageQC.csv"),quote=FALSE)
write.csv(qcmat,file=file.path(projectname,"qcmat.csv"),quote=FALSE,row.names=FALSE)
write.csv(data.frame(gc=gc_qc),file=file.path(projectname,"gc_qc.csv"),quote=FALSE,row.names=FALSE)
saveRDS(ref_qc,file=file.path(projectname,"ref_qc.rds"))
write.table(sampname_qc,file=file.path(projectname,"sampname_qc.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

Y.nonzero<-Y_qc[apply(Y_qc,1,function(x)!any(x==0)),,drop=FALSE]
if(!nrow(Y.nonzero)) stop("No exons without zero counts after QC; loosen QC thresholds.")
pseudo.sample<-apply(Y.nonzero,1,function(x)exp(mean(log(x))))
N<-apply(apply(Y.nonzero,2,function(x)x/pseudo.sample),2,median)
write.csv(data.frame(sample=names(N),N=as.numeric(N)),file=file.path(projectname,"library_size_factor.csv"),quote=FALSE,row.names=FALSE)

chrs<-as.character(unique(seqnames(ref_qc)))
chr_num<-suppressWarnings(as.integer(chrs))
chrs<-chrs[order(is.na(chr_num),chr_num,chrs)]
write.table(chrs,file=file.path(projectname,"chrs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

message("Prep done. Exons QC: ",nrow(Y_qc)," samples: ",ncol(Y_qc)," chrs: ",paste(chrs,collapse=","))
