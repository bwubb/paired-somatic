###facets-snakemake.R

library(facets)
library(argparse)

p <- ArgumentParser()
p$add_argument("--id", help="Tumor id.")
p$add_argument("--input", help="snp-pileup csv.gz matrix.")
p$add_argument("--cval", type="integer", default=150, help="procSample cval (coarser segmentation).")
p$add_argument("--ndepth", type="integer", default=25, help="Minimum normal depth in preProcSample.")
p$add_argument("--gbuild", default="hg38", help="Genome build for GC correction (hg19/hg38/...).")
args <- p$parse_args()

datafile <- args$input
cval <- args$cval
ndepth <- args$ndepth
gbuild <- args$gbuild
outpath <- dirname(normalizePath(datafile))

rcmat <- readSnpMatrix(datafile)
xx <- preProcSample(rcmat, ndepth=ndepth, gbuild=gbuild)
oo <- procSample(xx, cval=cval)
fit <- emcncf(oo)

pdf(file.path(outpath, "copynumber_profile.pdf"))
plotSample(x=oo, emfit=fit, sname=args$id)
dev.off()

pdf(file.path(outpath, "fit_diagnostic.pdf"))
logRlogORspider(oo$out, oo$dipLogR)
dev.off()

write.csv(
  data.frame(purity=fit$purity, ploidy=fit$ploidy, dipLogR=oo$dipLogR),
  file=file.path(outpath, "purity_ploidy.csv"),
  row.names=FALSE,
  quote=FALSE
)

flags <- oo$flags
if (is.null(flags) || length(flags) == 0) {
  flags <- "OK"
}
writeLines(as.character(flags), con=file.path(outpath, "flags.txt"))

write.csv(fit$cncf, file=file.path(outpath, "segmentation_cncf.csv"), row.names=FALSE, quote=FALSE)

segments.txt <- data.frame(
  chromosome=fit$cncf$chrom,
  start.pos=fit$cncf$start,
  end.pos=fit$cncf$end,
  CNt=fit$cncf$tcn.em,
  A=fit$cncf$tcn.em - fit$cncf$lcn.em,
  B=fit$cncf$lcn.em
)
write.table(
  segments.txt,
  file=file.path(outpath, paste0(args$id, "_segments.txt")),
  sep="\t",
  row.names=FALSE,
  quote=FALSE
)
