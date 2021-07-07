library(CODEX2)

#normal_bams<-read.table("normal_bams.list")
#tumor_bams<-read.table("tumor_bams.list")

#all_bams<-as.character(rbind(normal_bams,tumor_bams)$V1)

if (exists("snakemake")) {
    all_bams<-as.character(snakemake@input)
    projectname<-snakemake@params[['project']] #should be dir and not a prefix
    bedFile<-snakemake@params[['bed']]
    print(bedFile)
}#else
head(all_bams)
tail(all_bams)

sampname<-gsub('.ready.bam','',basename(all_bams))
#projectname is not accurate, its really codex2 working dir

bambedObj<-getbambed(bamdir=all_bams,bedFile=bedFile,sampname=sampname,projectname=projectname)
bamdir<-bambedObj$bamdir;sampname<-bambedObj$sampname
ref<-bambedObj$ref;projectname<-bambedObj$projectname

gc<-getgc(ref)
mapp<-getmapp(ref)
values(ref)<-cbind(values(ref),DataFrame(gc,mapp))

coverageObj<-getcoverage(bambedObj,mapqthres=20)
Y<-coverageObj$Y
write.csv(Y,file=paste(projectname,'coverage.csv',sep='/'),quote=FALSE)

qcObj<-qc(Y,sampname,ref,cov_thresh=c(20,4000),length_thresh=c(20,2000),mapp_thresh=0.9,gc_thresh=c(20,80))


Y_qc<-qcObj$Y_qc
write.csv(Y_qc,file=paste(projectname,'coverageQC.csv',sep='/'),quote=FALSE)

sampname_qc<-qcObj$sampname_qc
ref_qc<-qcObj$ref_qc
qcmat<-qcObj$qcmat
gc_qc<-ref_qc$gc
write.csv(gc_qc,file=paste(projectname,'gc_qc.csv',sep='/'),quote=FALSE)

#write.table(qcmat,file=paste(projectname,'qcmat.csv',sep='/'),sep='',quote=FALSE,row.names=FALSE)
write.csv(qcmat,file=paste(projectname,'qcmat.csv',sep='/'),quote=FALSE,row.names=FALSE)

Y.nonzero<-Y_qc[apply(Y_qc,1,function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N<-apply(apply(Y.nonzero,2,function(x){x/pseudo.sample}),2,median)

write.csv(N,file=paste(projectname,'library_size_factor.csv',sep='/'),quote=FALSE,row.names=FALSE)
