library(ASCAT)
library(argparse)

parser <- ArgumentParser(description='Run ASCAT')
parser$add_argument('--tumor-logr', required=TRUE, help='Tumor LogR file')
parser$add_argument('--tumor-baf', required=TRUE, help='Tumor BAF file')
parser$add_argument('--germline-logr', required=TRUE, help='Germline LogR file')
parser$add_argument('--germline-baf', required=TRUE, help='Germline BAF file')
parser$add_argument('--outdir', required=TRUE, help='Output directory')
parser$add_argument('--gc-file', required=TRUE, help='GC content file')
parser$add_argument('--rdata', required=TRUE, help='Output RData file')

#Note: They switched back to Tumor...
args <- parser$parse_args()

ascat.bc = ascat.loadData(
    Tumor_LogR_file=args$tumor_logr,
    Tumor_BAF_file=args$tumor_baf,
    Germline_LogR_file=args$germline_logr,
    Germline_BAF_file=args$germline_baf,
    gender='XX',
    genomeVersion="hg38",
    isTargetedSeq=T)

ascat.plotRawData(ascat.bc, img.dir=args$outdir, img.prefix="Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile=args$gc_file)
ascat.plotRawData(ascat.bc, img.dir=args$outdir, img.prefix="After_correction_")
ascat.bc = ascat.aspcf(ascat.bc, out.dir=args$outdir, penalty=25)
ascat.plotSegmentedData(ascat.bc, img.dir=args$outdir)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments=T, img.dir=args$outdir)
QC = ascat.metrics(ascat.bc, ascat.output)
write.table(QC, file=file.path(args$outdir, "QC.txt"), sep=" ")
save(ascat.bc, ascat.output, QC, file=args$rdata)
