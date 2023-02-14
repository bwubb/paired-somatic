###facets-snakemake.R

library("facets")


p<-ArgumentParser()
p$add_argument('--id',help='Tumor id.')
p$add_argument('--input',help='Input seqz file.')
p$add_argument('--cval',default=150,help="FACETS cval")
#args=p$parse_args(c("--output_dir","codex2","--bed_file","/home/bwubb/resources/Bed_files/MelanomaTargeted_v3.S3094091.Covered.bed","--bam_table","bam.table","--normal_list","normal.list"))
args<-p$parse_args()

#init
datafile=args$input[[1]]
cval=args$params[["cval"]]
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

#seq.cols.needed = c("chromosome", "start.pos", "end.pos", "CNt", "A", "B")
segments.txt=data.frame("chromosome"=fit$cncf$chrom, "start.pos"=fit$cncf$start,"end.pos"=fit$cncf$end,"CNt"=fit$cncf$tcn.em,"A"=fit$cncf$tcn.em-fit$cncf$lcn.em,"B"=fit$cncf$lcn.em)
write.table(segments.txt,file=paste(outpath,pasete0(args$id,"_segments.txt"),sep="/"),sep="\t",row.names=FALSE,quote=FALSE)
