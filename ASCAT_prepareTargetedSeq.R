library(ASCAT)

setClass(
  "MockSnakemake",
  representation(
    input = "list",
    params = "list",
    threads = "numeric"
  )
)

snakemake <- new("MockSnakemake",
  input = list(worksheet = "data/work/ASCAT/worksheet.tsv"),
  params = list(
    outdir="data/work/ASCAT",
    alleles = "/home/bwubb/projects/021-POSH_FFPE/G1000_allelesAll_hg38/G1000_alleles_hg38_chr",
    bed = "/home/bwubb/resources/Bed_files/SureSelect-Exon_v7.S31285117.GRCh38.bed"
  ),
  threads = 8  # Replace with the desired number of threads
)

ascat.prepareTargetedSeq(
       Worksheet=snakemake@input[["worksheet"]],
       Workdir=snakemake@params[["outdir"]],
       alleles.prefix=snakemake@params[["alleles"]],
       BED_file=snakemake@params[["bed"]],
       allelecounter_exe="/home/bwubb/software/alleleCount/bin/alleleCounter",
       genomeVersion="hg38",
       nthreads = snakemake@threads,
       minCounts = 10,
       is_chr_based = F,
       chrom_names = c(1:22, "X"),
       min_base_qual = 20,
       min_map_qual = 35,
       ref.fasta = NA,
       plotQC = T
     )
