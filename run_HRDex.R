#!/usr/bin/env Rscript
library(HRDex)
library(argparse)

p<-ArgumentParser()
p$add_argument('-i','--infile',help='Path to input file.')
p$add_argument('-o','--outfile',help="Path to output file.")
p$add_argument('--tumor',help='Tumor id.')
p$add_argument('--build',help='Genome build; {grch37,grch38}.')

args<-p$parse_args()

seq.dat <- read.table(args$infile, header = TRUE)
ref <- tolower(args$build)
sample.id <- args$tumor

#column names for facets
#csv text or csv

seq.dat <- preprocessHRD( seq.dat, ref )#needs required fields
#Does not check for A or B NA values, which facet may have.

CN.dat <- getCNt(seq.dat )

#outputting all values for now.
HRD.score <- getHRD.Score( seq.dat, CN.dat, scaleTotal = FALSE)
HRD.LST  <- getLST(  seq.dat )
HRD.LOH   <- getLOH(  seq.dat )
HRD.NTAIr <- getNTAI.raw( seq.dat )
HRD.NTAIm <- getNTAI.norm( seq.dat, CN.dat )

output=data.frame("TumorID"=args$tumor,"HRD.NTAIr"=HRD.NTAIr,"HRD.NTAIm"=HRD.NTAIm,"HRD.LOH"=HRD.LOH,"HRD.LST"=HRD.LST,"HRD.score"=HRD.score)
write.table(output,file=args$outfile,sep=",",row.names=FALSE,col.names=TRUE)
