

suppressMessages(library(IRanges))
suppressMessages(library(dplyr))
suppressMessages(library(Hmisc))
suppressMessages(library(stringr))
#suppressMessages(library(binom))

expSomaticVAF=function(phi,N,n)
{
  return((phi*n)/(phi*N+2*(1-phi)))
}

expGermlineVAF=function(phi,N,n)
{
  return((phi*n+(1-phi))/(phi*N+2*(1-phi)))
}

keep_words=function(text)
{
  keep=c("loh","gain","amp","del","loss","neutral")
  words <- strsplit(text, ";")[[1]]
  txt <- paste(words[words %in% keep], collapse = ";")
  return(txt)
}

SomaticIntersect=function(segment_data,somatic_data,cellularity)
{
  somatic_out=data.frame()
  for( i in c(1:22, "X") )
  {
    somatic.tmp=somatic_data[somatic_data$Chr==i,]
    chr.segments=segment_data[segment_data$chromosome==i,]
    #overlap of variants and segments
    ind=findOverlaps(IRanges(somatic.tmp$Start,somatic.tmp$Start),IRanges(chr.segments$start.pos,chr.segments$end.pos))
    somatic.tmp=somatic.tmp[ind@from,]
    somatic.tmp$N=as.integer(chr.segments$CNt[ind@to])
    somatic.tmp$n=as.integer(chr.segments$A[ind@to])
    #where is seqzLOH
    somatic.tmp$Sequenza.LOH=chr.segments$Sequenza.LOH[ind@to]
    somatic.tmp$expVAF=expSomaticVAF(cellularity,somatic.tmp$N,somatic.tmp$n)
    somatic_out=rbind(somatic_out,somatic.tmp)
  }
  return(somatic_out)
}

GermlineIntersect=function(segment_data,germline_data,cellularity)
{
  germline_out=data.frame()
  for( i in c(1:22, "X") )
  {
    germline.tmp=germline_data[germline_data$Chr==i,]
    chr.segments=segment_data[segment_data$chromosome==i,]
    ind2=findOverlaps(IRanges(germline.tmp$Start,germline.tmp$Start),IRanges(chr.segments$start.pos,chr.segments$end.pos))
    germline.tmp=germline.tmp[ind2@from,]
    germline.tmp$N=as.integer(chr.segments$CNt[ind2@to])
    germline.tmp$n=as.integer(chr.segments$A[ind2@to])
    germline.tmp$Sequenza.LOH=chr.segments$Sequenza.LOH[ind2@to]
    germline.tmp$expVAF=expGermlineVAF(cellularity,germline.tmp$N,germline.tmp$n)
    germline_out=rbind(germline_out,germline.tmp)
  }
  return(germline_out)
}

read.annotated_file=function(filename)
{
  file.data=read.table(file,header=F,sep="\t",quote="")
  file.data=file.data[,c("TumorID","Chr","Start","Ref","Alt","GenomicRegion.refGene","ExonicFunc.refGene","Gene.refGene","Exon.refGene","NTChange.refGene","AAChange.refGene","Tumor_Zyg","Tumor_Total_Depth","Tumor_ALT_AlleleDepth","Tumor_ALT_AlleleFrac","CLNSIG")]
}

###If snakemake
suppressPackageStartupMessages(library("argparse"))
# create parser object
parser=ArgumentParser()
#Switch to sequenza output directly
#Would have to make intervals from it.
#Or add sequenza to bed into sequenza.snake
parser$add_argument('--segments_bed',type="character",help='')
#Switch to vcf in future.
parser$add_argument('--germline_csv',type="character",help='')
parser$add_argument('--somatic_csv',type="character",help='')
#
parser$add_argument('--cellularity_file',type="character",help='Sequenza, cellularity file. Dont use if passing cellularity directly')
parser$add_argument('--cellularity',type="double",default=1.0,help='Ceullarity. 0.0-1.0')

argv=parser$parse_args()
segments_bed=argv$segments_bed
somatic_csv=argv$somatic_csv
germline_csv=argv$germline_csv

if (file.exists(argv$cellularity_file))
{
    cellularity=as.numeric(read.table(file=argv$cellularity_file,sep="\t",header=TRUE)[2,1])
}else
{
    cellularity=argv$cellularity
}

print(cellularity)
chrOrder <-c((1:22),"X","Y","M")
#Load Segments
print(paste("Reading",segments_bed))
segments.data=read.table(segments_bed,header=F,sep="\t",quote="")
colnames(segments.data)=c("chromosome","start.pos","end.pos","name","score","strand")

segments.data$CNt=segments.data$score/100
segments.data$A=as.integer(str_match(segments.data$name,"A([0-9])+")[,2])
segments.data$seqz=unlist(lapply(as.character(segments.data$name),keep_words))
segments.data$Sequenza.LOH[!grepl("loh",segments.data$seqz)]="NO"
segments.data$Sequenza.LOH[grepl("loh",segments.data$seqz)]="YES"
#Where is allelic imbalance?
#segments.data$B=segments.data$CNt-segments.data$A
#AI=segments.data$B>0&segments.data$A>segments.data$B


#Load Variant annotation w/ Freq
print(paste("Reading",somatic_csv))
somatic.data=read.csv(file=somatic_csv,header=T,quote="")
somatic.data=somatic.data[,c("TumorID","Chr","Start","Ref","Alt","GenomicRegion.refGene","ExonicFunc.refGene","Gene.refGene","Exon.refGene","NTChange.refGene","AAChange.refGene","Tumor_Zyg","Tumor_Total_Depth","Tumor_ALT_AlleleDepth","Tumor_ALT_AlleleFrac","CLNSIG")]
somatic.data$Type="SOMATIC"
#maybe add NumPass,mutect2, etc.

print(paste("Reading",germline_csv))
germline.data=read.csv(file=germline_csv,header=T,quote="")
germline.data=germline.data[,c("TumorID","Chr","Start","Ref","Alt","GenomicRegion.refGene","ExonicFunc.refGene","Gene.refGene","Exon.refGene","NTChange.refGene","AAChange.refGene","Tumor_Zyg","Tumor_Total_Depth","Tumor_ALT_AlleleDepth","Tumor_ALT_AlleleFrac","CLNSIG")]
germline.data$Type="GERMLINE"
#Where are the HOM_REF?
#FILTER depth
#No HOM_REF in Normal

#####
somatic.intersect=SomaticIntersect(segments.data,somatic.data,cellularity)

somatic_conf=as.data.frame(binconf( sum(somatic.intersect$N==somatic.intersect$n),dim(somatic.intersect)[1]))
somatic.intersect$Variant.LOH="nonLOH" #GenomicFunc?
somatic.loh.i=somatic.intersect$Sequenza.LOH=="YES"
somatic.intersect[somatic.loh.i,][which(somatic.intersect[somatic.loh.i,]$Tumor_ALT_AlleleFrac>somatic_conf$Lower),"Variant.LOH"]="LOH"

#somatic.ai.i=(somatic.intersect$N-somatic.intersect$n)>0&somatic.intersect$n>(somatic.intersect$N-somatic.intersect$n)
#somatic.intersect[somatic.ai.i,][which(between(somatic.intersect[somatic.ai.i,]$Tumor_ALT_AlleleFrac,somatic_conf$Lower,somatic_conf$Upper)),"Variant.LOH"]="LOH"

#make better variable names


#####

germline.intersect=GermlineIntersect(segments.data,germline.data,cellularity)


germline_conf=as.data.frame(binconf( sum(germline.intersect$N==germline.intersect$n),dim(germline.intersect)[1]))
germline.intersect$Variant.LOH="nonLOH"
germline.loh.i=germline.intersect$Sequenza.LOH=="YES"
germline.intersect[germline.loh.i,][which(germline.intersect[germline.loh.i,]$Tumor_ALT_AlleleFrac>germline_conf$Lower),"Variant.LOH"]="LOH"

germline.intersect[germline.loh.i,][which(1-germline.intersect[germline.loh.i,]$Tumor_ALT_AlleleFrac>germline_conf$Lower),"Variant.LOH"]="LOH"

out.intersect=rbind(somatic.intersect,germline.intersect)%>%arrange(factor(Chr,chrOrder,ordered = TRUE),Start)


write.csv(out.intersect,file=paste(dirname(germline_csv),"variant_LOH.csv",sep="/"),quote = F,row.names = F)


#Pull out LOH only.
#Redo seperately allelic imbalance
#seqzLOH needs to be loh;amp/del/etc
#hidden, not written, columns. seqzLOH.b

#### 2nd pass
LOHData=function(segments.data,out.intersect)
{
  loh_count=data.frame()
  for( i in c(1:22, "X") )
  {
    test.tmp=out.intersect[out.intersect$Chr==i,]
    chr.segments=segments.data[segments.data$chromosome==i,]
    for(j in c(1:nrow(chr.segments)))
    {
      ind=findOverlaps(IRanges(chr.segments$start.pos[j],chr.segments$end.pos[j]),IRanges(test.tmp$Start,test.tmp$Start))
      tmp=test.tmp[ind@to,]
      loh_count=rbind(loh_count,data.frame(Chr=i,
                                           Start=chr.segments$start.pos[j],
                                           End=chr.segments$end.pos[j],
                                           CNV.type=chr.segments$seqz[j],
                                           Length=chr.segments$end.pos[j]-chr.segments$start.pos[j],
                                           CN.t=chr.segments$CNt[j],
                                           CN.a=chr.segments$A[j],
                                           Variant.Total=nrow(tmp),
                                           LOH.count=sum(tmp$Variant.LOH=="LOH"),
                                           PCT=sum(tmp$Variant.LOH=="LOH")/nrow(tmp),
                                           LOH.SOMATIC=sum(tmp$Variant.LOH=="LOH"&tmp$Type=="SOMATIC"),
                                           NO.LOH.SOMATIC=sum(tmp$Variant.LOH!="LOH" & tmp$Type=="SOMATIC"),
                                           LOH.GERMLINE=sum(tmp$Variant.LOH=="LOH" & tmp$Type=="GERMLINE"),
                                           NO.LOH.GERMLINE=sum(tmp$Variant.LOH!="LOH" & tmp$Type=="GERMLINE")))
      #Good place to plot variants, Color by SOMATIC/GERMLINE. Add a few specific names of important mutations
    }
  }
  return(loh_count)
}

loh_count=LOHData(segments.data, out.intersect)

write.csv(loh_count,file=paste(dirname(germline_csv),"segments_LOH.csv",sep="/"),quote = F,row.names = F)
#How well do the loh calls fit? I can weight segments by size. Get % of total
#number agree somatic/germline, number disagree somatic/germline
