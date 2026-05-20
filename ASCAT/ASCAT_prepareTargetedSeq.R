library(ASCAT)
library(argparse)

parser <- ArgumentParser(description='Prepare targeted sequencing data for ASCAT')
parser$add_argument('--worksheet', required=TRUE, help='Path to sample worksheet TSV')
parser$add_argument('--outdir', required=TRUE, help='Output directory')
parser$add_argument('--alleles-prefix', required=TRUE, help='Prefix path for alleles files')
parser$add_argument('--bed', required=TRUE, help='Path to target BED file')
parser$add_argument('--allelecounter', required=TRUE, help='Path to alleleCounter executable')
parser$add_argument('--threads', type='integer', default=8, help='Number of threads')
parser$add_argument('--min-counts', type='integer', default=10, help='Minimum read counts')
parser$add_argument('--min-base-qual', type='integer', default=20, help='Minimum base quality')
parser$add_argument('--min-map-qual', type='integer', default=35, help='Minimum mapping quality')
#parser$add_argument('--ref-fasta', help='Reference FASTA file path', default=NA)

args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive=TRUE)
}

# Create alleleData directory structure
allele_data_dir <- file.path(args$outdir, "alleleData")
if (!dir.exists(allele_data_dir)) {
    dir.create(allele_data_dir, recursive=TRUE)
}
if (!dir.exists(file.path(allele_data_dir, "Cleaned"))) {
    dir.create(file.path(allele_data_dir, "Cleaned"), recursive=TRUE)
}

ascat.prepareTargetedSeq(
    Worksheet=args$worksheet,
    Workdir=args$outdir,
    alleles.prefix=args$alleles_prefix,
    BED_file=args$bed,
    allelecounter_exe=args$allelecounter,
    genomeVersion="hg38",
    nthreads=args$threads,
    minCounts=args$min_counts,
    is_chr_based=FALSE,#Need a way to handle this.
    chrom_names=c(1:22, "X"),
    min_base_qual=args$min_base_qual,
    min_map_qual=args$min_map_qual,
    plotQC=TRUE
)
