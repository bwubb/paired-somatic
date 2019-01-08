
library(CODEX2)

#args = commandArgs(trailingOnly=TRUE)

#normal_bams<-read.table("normal_bams.list")
#tumor_bams<-read.table("tumor_bams.list")
normal_bams=snakemake@input[['normals']]
tumor_bams=snakemake@input[['tumors']]

all_bams<-as.character(rbind(normal_bams,tumor_bams)$V1)
head(all_bams)
tail(all_Bams)

sampname<-gsub(".ready.bam","",basename(all_bams))
bedFile<-snakemake@params[['bed']]

#Replace bamdir with all bams
bambedObj<-getbambed(bamdir=all_bams,bedFile=bedFile,sampname=sampname,projectname="tumor_normal")
bamdir<-bambedObj$bamdir
sampname<-bambedObj$sampname
ref<-bambedObj$ref
projectname<-bambedObj$projectname

#OR SNAKEMAKE
#Have to provide a list of tumor files and a list of normal files
#Follow Yuchao's bullshit with the combined list being bamFile
#get index of normals
#bambedObj <- getbambed(bamdir = snakemake@params[[]], bedFile = bedFile, sampname = sampname, projectname = '')


##Getting GC content and mappability
#Obtain GC contents for each exon/target/window.
#gc<-getgc(ref)
#mapp<-getmapp(ref)
#values(ref)<-cbind(values(ref),DataFrame("gc"=gc,"mapp"=mapp))

#coverageObj <- getcoverage(bambedObj, mapqthres = 20)
#Y <- coverageObj$Y
#write.csv(Y, file = paste(projectname, '_coverage.csv', sep=''), quote = FALSE)

##Quality control
#qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 4000),length_thresh = c(20, 2000), mapp_thresh = 0.9,gc_thresh = c(20, 80))
#Y_qc<-qcObj$Y_qc
#write.csv(Y_qc,file=paste(projectname,'_coverageQC.csv',sep=''),quote=FALSE)

#Estimate library size factor based on genome-wide read depth after QC.
#Y.nonzero<-Y_qc[apply(Y_qc,1,function(x){!any(x==0)}),]
#pseudo.sample<-apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
#N<-apply(apply(Y.nonzero,2,function(x){x/pseudo.sample}),2,median)
#write.csv(N,file=paste(projectname,'_N','.csv',sep=''),quote=FALSE,row.names=FALSE)

#sampname_qc<-qcObj$sampname_qc
#ref_qc<-qcObj$ref_qc
#qcmat<-qcObj$qcmat
#gc_qc<-ref_qc$gc
#write.csv(gc_qc,file=paste(projectname,'_gc_qc.csv',sep=''),quote=FALSE)
#write.table(qcmat,file=paste(projectname,'_qcmat','.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

### Running CODEX2 with Negative Control Samples###

Y_qc <- read.csv(paste(projectname,'codex2_coverageQC.csv',sep='/'),row.names=1,check.names=FALSE)
Y_qc<-as.matrix(Y_qc)
message("nrow of Y_qc ", nrow(Y_qc))

qcmat<-read.table(file=paste(projectname,'codex2_qcmat.',sep='/'),sep='\t',header=TRUE)
ind<-qcmat$pass==TRUE
ref_qc<-ref[ind]
message('Size of ref_qc ', length(ref_qc))
#names?

gc<-getgc(ref_qc)
mapp<-getmapp(ref_qc)
values(ref_qc)<-cbind(values(ref_qc), DataFrame(gc, mapp))
gc_qc<-read.csv(paste(projectname,'codex2_gc_qc.csv',sep='/'),row.names=1,check.names=FALSE)
gc_qc<-as.numeric(gc_qc[,'x'])
message('Length of gc_qc ', length(gc_qc))

N<-read.csv(paste(projectname,'library_size_factor.csv',sep='/'))
N<-as.numeric(N[,'x'])

###break###
norm_index<-which(colnames(Y_qc) %in% gsub(".ready.bam","",basename(as.character(normal_bams$V1))))

sampname_qc<-colnames(Y_qc)

#chr<-args[1]
chr<-snakemake.params[["chr"]]
message("chr",chr)

chr.index<-which(seqnames(ref_qc)==chr)
message('Length of chr.index ', length(chr.index))

message('Length of gc_qc[chr.index] ,' gc_qc[chr.index])

message('nrow of Y_qc[chr.index,] ', nrow(Y_qc[chr.index,]))

tail(Y_qc[chr.index,])

normObj<-normalize_codex2_ns(Y_qc=Y_qc[chr.index,],gc_qc=gc_qc[chr.index],K=1:10,norm_index=norm_index,N=N)

Yhat.ns<-normObj$Yhat
Yhat.ns
fGC.hat.ns<-normObj$fGC.hat
beta.hat.ns<-normObj$beta.hat
g.hat.ns<-normObj$g.hat
h.hat.ns<-normObj$h.hat

AIC.ns<-normObj$AIC
BIC.ns<-normObj$BIC
RSS.ns<-normObj$RSS

#Choose the number of latent Poisson factors. BIC is used as the model selection metric by default.
choiceofK(AIC.ns,BIC.ns,RSS.ns,K=1:10,filename=paste(projectname,'/codex2_chr',chr,'_ns_choiceofK.pdf',sep=''))


####break if you can recover Yhat.ns and BIC.ns
### Running segmentation by CODEX2 ###
BIC.ns
optK=which.max(BIC.ns)
message('optK ',optK)

#In segmentation step, use "fraction" mode for somatic CNA detection (cancer is heterogenous)
finalcall.CBS<-segmentCBS(Y_qc[chr.index,],Yhat.ns,optK=which.max(BIC.ns),K=1:10,sampname_qc=sampname_qc,ref_qc=ranges(ref_qc),chr=chr,lmax=400,mode="fraction")
write.table(finalcall.CBS,file=paste(projectname,'/codex2_chr',chr,'.K_',optK,'.segments.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

#Post-segmentation pruning and filtering are recommended based on CNV length (filter1), length per exon (filter2), likelihood ratio (filter3), and number of exons (filter4).

filter1<-finalcall.CBS$length_kb<=200
filter2<-finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
finalcall.CBS.filter<-finalcall.CBS[filter1&filter2,]

filter3<-finalcall.CBS.filter$lratio>40
filter4<-(finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
finalcall.CBS.filter=finalcall.CBS.filter[filter3|filter4,]

write.table(finalcall.CBS.filter,file=projectname,'/codex2_chr',chr,'.K_',optK,'.segments.filtered.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
