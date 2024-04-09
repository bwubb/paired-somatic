library(ASCAT)
#Note: The ASCAT package uses 'tumour' whereas I am American.

ascat.prepareHTS(
    tumourseqfile=snakemake@input[["tumor"]],
    normalseqfile=snakemake@input[["normal"]],
    tumourname=snakemake@params[["tumor"]],
    normalname=snakemake@params[["normal"]],
    allelecounter_exe="/home/bwubb/software/alleleCount/bin/alleleCounter",
    gender="XX",
    genomeVersion="hg38",
    BED_file=snakemake@params[["bed"]],
    alleles.prefix=snakemake@params[["alleles"]],
    loci.prefix=snakemake@params[["loci"]],
    nthreads=snakemake@threads,
    tumourLogR_file=snakemake@output[["tumorLogR"]],
    tumourBAF_file=snakemake@output[["tumorBAF"]],
    normalLogR_file=snakemake@output[["germlineLogR"]],
    normalBAF_file=snakemake@output[["germlineBAF"]]
    )

#mv allele fequencies files
