############################
#Brad Wubbenhorst
#bwubb@pennmedicine.upenn.edu
#Dec. 2018

library(sequenza)
library(argparse)
#data.file <- system.file("data", "example.seqz.txt.gz", package = "sequenza")

p<-ArgumentParser()
p$add_argument('--id',help='Tumor id.')
p$add_argument('--input',help='Input seqz file.')
p$add_argument('--outdir',help="Path to put all output files.")
p$add_argument('--threads',default=1,help="Number of cores")
#args=p$parse_args(c("--output_dir","codex2","--bed_file","/home/bwubb/resources/Bed_files/MelanomaTargeted_v3.S3094091.Covered.bed","--bam_table","bam.table","--normal_list","normal.list"))
args<-p$parse_args()

tumor<-sequenza.extract(args$input)
CP<-sequenza.fit(tumor,mc.cores=args$threads)
sequenza.results(sequenza.extract=tumor,cp.table=CP,sample.id=args$id,out.dir=args$outdir)
