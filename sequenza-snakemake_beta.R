############################
#Brad Wubbenhorst
#bwubb@pennmedicine.upenn.edu
#Dec. 2018

library(sequenza)
# predefined data about chromosome size, centromere, and telomere location
# chromosome size, centromere and telomere locations (in hg19/GRCh37)
ref.dat = data.frame( chromosome = c(seq(1:22), "X", "Y"),
      centromere.start = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166,
                           58054331, 43838887, 47367679, 39254935, 51644205, 34856694, 
                           16000000, 16000000, 17000000, 35335801, 22263006, 15460898,
                           24681782, 26369569, 11288129, 13000000, 58632012, 10104553),
      centromere.end = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166,
                         61054331, 46838887, 50367679, 42254935, 54644205, 37856694,
                         19000000, 19000000, 20000000, 38335801, 25263006, 18460898,
                         27681782, 29369569, 14288129, 16000000, 61632012, 13104553),
      p.telomere.end = rep(10000, 24),
      q.telomere.start = c(249240621, 243189373, 198012430, 191144276, 180905260,
                           171105067, 159128663, 146354022, 141203431, 135524747,
                           134996516, 133841895, 115159878, 107339540, 102521392,
                           90344753, 81185210, 78067248, 59118983, 63015520,
                           48119895, 51294566, 155260560, 59363566),
      chr.size = c(249250621, 243199373,  198022430, 191154276, 180915260,
                   171115067, 159138663, 146364022, 141213431, 135534747,
                   135006516, 133851895, 115169878, 107349540, 102531392,
                   90354753, 81195210, 78077248, 59128983, 63025520,
                    48129895, 51304566, 155270560, 59373566))

#input="seqz.small.gz"
robustExtract<-function(x) {tryCatch(sequenza.extract(input,chromosome.list=x),error=function(e) { print(paste("Extract Failed: Chromosome", x)); NaN})}

#results<-vector("list",23)
#names(results)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
results<-vector("list",23)
names(results)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")

for (chr in names(results))
    {
    results[chr]<-robustExtract(chr)
    #Warning: number of items to replace is not a multiple of replacement length
    }

if (anyNA(results))
    {
    print(is.na(results))
    break.list<-vector("list",3)
    names(break.list)<-c("chrom","start.pos","end.pos")
    for (chr in names(results))
    {
        if (is.na(results[chr]))
        {
            #add p arm and q arm from ref.data
            break.list["chrom"]<-c(break.list["chrom"],chr)
            break.list["star.pos"]<-c(break.list["star.pos"],"1")
            break.list["end.pos"]
    }

#if yes take FALSE  results, take the first 3 columns: chromosome, start, end and save in a object/file (name the columns chrom, start.pos, end.pos)
#open the file/object and manually add a line for chromosome 17 (eg from start to end of the chromosome) 
#maybe loop? for chr in ... if results[chr]==NaN then make_segments(chr)

  #-R "select[hname!=node047.hpc.local]"
tumor<-sequenza.extract(snakemake@input[[1]])
CP<-sequenza.fit(tumor,mc.cores=snakemake@threads)
sequenza.results(sequenza.extract=tumor,cp.table=CP,sample.id=snakemake@wildcards[["tumor"]],out.dir=snakemake@params[["outdir"]])