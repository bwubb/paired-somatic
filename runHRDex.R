#!/usr/bin/env Rscript
library(HRDex)
library(argparse)

p <- ArgumentParser()
p$add_argument('-i', '--infile', help = 'Path to input BED file.')
p$add_argument('-o', '--outfile', help = "Path to output file.")
p$add_argument('--tumor', help = 'Tumor id.')
p$add_argument('--build', help = 'Genome build; {grch37, grch38, hg19, hg38}.')

args <- p$parse_args()

#test data
args<-list(infile="05217-002-DZ1A.call.bed",build="grch38",tumor="05217-002-DZ1A")

# Mapping hg19 to grch37 and hg38 to grch38
genome_build_map <- list(grch37 = "grch37", grch38 = "grch38", hg19 = "grch37", hg38 = "grch38")
ref <- tolower(genome_build_map[[tolower(args$build)]])

if (is.null(ref)) {
  stop("Invalid genome build. Please use one of: grch37, grch38, hg19, hg38")
}

# Reading the BED file and specifying column names
seq.dat <- read.table(args$infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(seq.dat) <- c("chromosome", "start.pos", "end.pos", "name", "score", "strand")

# Extracting A and B values from the 'name' column and calculating CNt
seq.dat$A <- as.numeric(sub(".*;A(\\d+);.*", "\\1", seq.dat$name))
seq.dat$B <- as.numeric(sub(".*;B(\\d+).*", "\\1", seq.dat$name))

seq.dat <- na.omit(seq.dat)
seq.dat <- seq.dat[!is.na(seq.dat$A) & !is.na(seq.dat$B),]
if (nrow(seq.dat) > 0) {
    seq.dat$CNt <- seq.dat$A + seq.dat$B
  } else {
    warning("No valid rows found for CNt calculation in file: ", file_path)
    return(NULL) # Skip further processing for this file
  }

# Ensuring chromosome format is consistent for analysis
seq.dat$chromosome <- as.character(seq.dat$chromosome)
seq.dat$chromosome <- ifelse(seq.dat$chromosome == "23", "X", ifelse(seq.dat$chromosome == "24", "Y", seq.dat$chromosome))

sample.id <- args$tumor

# Now seq.dat includes the required fields:
seq.dat <- preprocessHRD(seq.dat, ref) # Make sure preprocessHRD is ready to accept this structure

CN.dat <- getCNt(seq.dat)

# Computing HRD scores
HRD.score <- getHRD.Score(seq.dat, CN.dat, scaleTotal = FALSE)
HRD.LST <- getLST(seq.dat)
HRD.LOH <- getLOH(seq.dat)
HRD.NTAIr <- getNTAI.raw(seq.dat)
HRD.NTAIm <- getNTAI.norm(seq.dat, CN.dat)

# Preparing output
output <- data.frame(TumorID = args$tumor, HRD.NTAIr = HRD.NTAIr, HRD.NTAIm = HRD.NTAIm, HRD.LOH = HRD.LOH, HRD.LST = HRD.LST, HRD.score = HRD.score)
write.csv(output, file = args$outfile, row.names = FALSE, col.names = TRUE)
