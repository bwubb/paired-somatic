
library(CODEX2)
args = commandArgs(trailingOnly=TRUE)

#normal_samples<-read.table(file=snakemake@params[['normals']])
#projectname<-snakemake@params[['project']]
normal_samples<-read.table(file='normals.S07604715.list')
projectname<-'data/work/BasserExome_S07604715/codex2'
chr<-args[1]
bedFile<-'/home/bwubb/resources/Bed_files/SureSelect-Exon_v6+COSMIC.S07604715.Covered.bed'
exomtarg <- read.table(bedFile, sep = "\t")
ref <- GRanges(seqnames=exomtarg[,1],ranges=IRanges(start=exomtarg[,2], end=exomtarg[,3]))
ref <- sort(ref)

### Running CODEX2 with Negative Control Samples###

Y_qc <- read.csv(paste(projectname,'coverageQC.csv',sep='/'),row.names=1,check.names=FALSE)
Y_qc<-as.matrix(Y_qc)
message("nrow of Y_qc ", nrow(Y_qc))

qcmat<-read.csv(file=paste(projectname,'qcmat.csv',sep='/'),header=TRUE)
#restrict chr here at qcmat
ind<-qcmat$pass==TRUE
ref_qc<-ref[ind]
message('Size of ref_qc ', length(ref_qc))
#names?

#I need to restrict to chr
gc<-getgc(ref_qc)
mapp<-getmapp(ref_qc)
values(ref_qc)<-cbind(values(ref_qc), DataFrame(gc, mapp))
gc_qc<-read.csv(paste(projectname,'gc_qc.csv',sep='/'),row.names=1,check.names=FALSE)
gc_qc<-as.numeric(gc_qc[,'x'])
message('Length of gc_qc ', length(gc_qc))

N<-read.csv(paste(projectname,'library_size_factor.csv',sep='/'))
N<-as.numeric(N[,'x'])

###break###
norm_index<-which(colnames(Y_qc) %in% normal_samples$V1)

sampname_qc<-colnames(Y_qc)

#chr<-args[1]
#chr<-snakemake.params[["chr"]]
message("chr",chr)

#If chr is pre-subset then this subsetting will change
chr.index<-which(seqnames(ref_qc)==chr)
message('Length of chr.index ', length(chr.index))

message('Length of gc_qc[chr.index] ', length(gc_qc[chr.index]))

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
message(paste('Opt K',optK,sep=' '))
write.table(finalcall.CBS,file=paste(projectname,'/chr',chr,'.codex2.segments.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

#Post-segmentation pruning and filtering are recommended based on CNV length (filter1), length per exon (filter2), likelihood ratio (filter3), and number of exons (filter4).

filter1<-finalcall.CBS$length_kb<=200
filter2<-finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
finalcall.CBS.filter<-finalcall.CBS[filter1&filter2,]

filter3<-finalcall.CBS.filter$lratio>40
filter4<-(finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
finalcall.CBS.filter=finalcall.CBS.filter[filter3|filter4,]

write.table(finalcall.CBS.filter,file=paste(projectname,'/chr',chr,'.codex2.segments.filtered.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
