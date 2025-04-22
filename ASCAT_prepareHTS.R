library(ASCAT)
library(argparse)

parser <- ArgumentParser(description='Prepare HTS data for ASCAT')
parser$add_argument('--tumor-bam', required=TRUE, help='Tumor BAM file')
parser$add_argument('--normal-bam', required=TRUE, help='Normal BAM file')
parser$add_argument('--tumor-name', required=TRUE, help='Tumor sample name')
parser$add_argument('--normal-name', required=TRUE, help='Normal sample name')
parser$add_argument('--outdir', required=TRUE, help='Output directory')
parser$add_argument('--bed', required=TRUE, help='Target BED file')
parser$add_argument('--alleles-prefix', required=TRUE, help='Prefix for alleles files')
parser$add_argument('--loci-prefix', required=TRUE, help='Prefix for loci files')
parser$add_argument('--allelecounter', required=TRUE, help='Path to alleleCounter executable')
parser$add_argument('--threads', type='integer', default=1, help='Number of threads')

args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive=TRUE)
}

# Define output files using file.path
tumor_logr <- file.path(args$outdir, "Tumor_LogR.txt")
tumor_baf <- file.path(args$outdir, "Tumor_BAF.txt")
normal_logr <- file.path(args$outdir, "Germline_LogR.txt")
normal_baf <- file.path(args$outdir, "Germline_BAF.txt")

#Note: The ASCAT package uses 'tumour' whereas I am American.
ascat.prepareHTS(
    tumourseqfile=args$tumor_bam,
    normalseqfile=args$normal_bam,
    tumourname=args$tumor_name,
    normalname=args$normal_name,
    allelecounter_exe=args$allelecounter,
    gender="XX",
    genomeVersion="hg38",
    BED_file=args$bed,
    alleles.prefix=args$alleles_prefix,
    loci.prefix=args$loci_prefix,
    nthreads=args$threads,
    tumourLogR_file=tumor_logr,
    tumourBAF_file=tumor_baf,
    normalLogR_file=normal_logr,
    normalBAF_file=normal_baf
)

#mv allele fequencies files
