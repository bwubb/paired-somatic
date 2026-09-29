#!/usr/bin/env Rscript
#Germline CODEX2 prep.
#  GRCh38 -> BSgenome.Hsapiens.NCBI.GRCh38 (no chr)
#  hg38   -> BSgenome.Hsapiens.UCSC.hg38 (chr*)
#
#Rscript codex2_prep.R \
#  --bam-table bam.table \
#  --samples data/work/codex2/samples.list \
#  --bed /path/to/Covered.bed \
#  --genome GRCh38 \
#  [--outdir data/work/codex2] \
#  [--mapq 20] \
#  [--mappability none|getmapp]

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
if(is.null(args[["bam-table"]])||is.null(args$samples)||is.null(args$bed)||is.null(args$genome)){
    stop("Required: --bam-table --samples --bed --genome GRCh38|hg38")
}
bam_table<-args[["bam-table"]]
sample_list<-args$samples
bedFile<-args$bed
outdir<-if(!is.null(args$outdir)) args$outdir else "data/work/codex2"
genome_key<-args$genome
mapq<-if(!is.null(args$mapq)) as.integer(args$mapq) else 20L
mapp_mode<-if(!is.null(args$mappability)) args$mappability else "none"
if(!mapp_mode%in%c("none","getmapp")) stop("--mappability must be none or getmapp (got: ",mapp_mode,")")

suppressPackageStartupMessages(library(CODEX2))

if(identical(genome_key,"GRCh38")){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
    genome<-BSgenome.Hsapiens.NCBI.GRCh38
    expect_chr_prefix<-FALSE
}else if(identical(genome_key,"hg38")){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    genome<-BSgenome.Hsapiens.UCSC.hg38
    expect_chr_prefix<-TRUE
}else{
    stop("--genome must be GRCh38 or hg38 (got: ",genome_key,")")
}

if(!file.exists(bam_table)) stop("bam.table not found: ",bam_table)
if(!file.exists(sample_list)) stop("sample list not found: ",sample_list)
if(!file.exists(bedFile)) stop("bed not found: ",bedFile)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

bams<-read.delim(bam_table,header=FALSE,sep="\t",stringsAsFactors=FALSE,col.names=c("sample","bam"))
bampath<-setNames(bams$bam,bams$sample)

samples<-scan(sample_list,what=character(),quiet=TRUE)
samples<-samples[nzchar(samples)&!grepl("^#",samples)]
if(!length(samples)) stop("empty sample list")
missing_id<-setdiff(samples,names(bampath))
if(length(missing_id)) stop("samples not in bam.table:\n",paste(missing_id,collapse="\n"))

bamdir<-unname(bampath[samples])
sampname<-samples
missing_bam<-bamdir[!file.exists(bamdir)]
if(length(missing_bam)) stop("missing BAM(s):\n",paste(missing_bam,collapse="\n"))
if(anyDuplicated(sampname)) stop("duplicate sample IDs")

message("Samples: ",length(sampname))
message("BED: ",bedFile)
message("Outdir: ",outdir)
message("genome: ",genome_key)
message("mappability: ",mapp_mode)

bambedObj<-getbambed(bamdir=bamdir,bedFile=bedFile,sampname=sampname,projectname=outdir)
bamdir<-bambedObj$bamdir
sampname<-bambedObj$sampname
ref<-bambedObj$ref
projectname<-bambedObj$projectname

seqlev<-as.character(seqlevels(ref))
has_chr<-any(grepl("^chr",seqlev))
if(expect_chr_prefix&&!has_chr){
    stop("genome=hg38 expects chr* contigs in BED, but BED has none. Check resources.targets_bed.")
}
if(!expect_chr_prefix&&has_chr){
    stop("genome=GRCh38 expects no chr prefix, but BED has chr* (e.g. ",grep("^chr",seqlev,value=TRUE)[1],").")
}

#CODEX2::getgc() assumes UCSC chr* names (looks up "chr1" etc). With GRCh38
#BED/BAM contigs ("1","2",...) that breaks against NCBI BSgenome. Compute GC
#ourselves with getSeq so seqnames must match the genome package as-is.
get_gc_pct<-function(ref,genome){
    chrs<-unique(as.character(seqnames(ref)))
    missing<-setdiff(chrs,seqnames(genome))
    if(length(missing)){
        stop("contigs in BED not in BSgenome (",genome_key,"): ",paste(utils::head(missing,10),collapse=","))
    }
    seqs<-Biostrings::getSeq(genome,ref)
    as.numeric(Biostrings::letterFrequency(seqs,letters="CG",as.prob=TRUE)*100)
}

message("Computing GC (direct getSeq; not CODEX2::getgc)...")
gc<-get_gc_pct(ref,genome)
if(identical(mapp_mode,"getmapp")){
    message("Computing mappability via getmapp()...")
    #getmapp has the same chr* assumption; only reliable for genome=hg38.
    if(!expect_chr_prefix){
        warning("getmapp is UCSC/chr*-oriented; with GRCh38 no-chr contigs falling back to mappability=none")
        mapp<-rep(1,length(ref))
    }else{
        mapp<-tryCatch(getmapp(ref,genome=genome),error=function(e){
            warning("getmapp failed (",conditionMessage(e),"); falling back to mappability=none")
            rep(1,length(ref))
        })
    }
}else{
    message("mappability=none: all targets set to 1 (no mappability filter).")
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
chr_num<-suppressWarnings(as.integer(sub("^chr","",chrs)))
chrs<-chrs[order(is.na(chr_num),chr_num,chrs)]
write.table(chrs,file=file.path(projectname,"chrs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)

message("Prep done. Exons QC: ",nrow(Y_qc)," samples: ",ncol(Y_qc)," chrs: ",paste(chrs,collapse=","))
