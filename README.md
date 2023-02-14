# paired-somatic
Collection of code for a paired somatic variant calling pipeline.

***

## Alginment
#### input: unaligned fastq files or ubam.
#### output: reference aligned bam per sample.

### BWA
   [https://github.com/lh3/bwa]  
   
   * bwa-mem alignment to refernce fasta

### samtools
  [https://github.com/samtools/samtools]  
  
  * samtools addreplacerg > fixmate > sort  
  * samtools markdup

### Metrics
   [https://broadinstitute.github.io/picard/command-line-overview.html]  
   
   * picardtools AlignmentMetrics and CollectHsMetrics
   * samtools markdup output

***

## Variant Discovery
#### input: reference aligned bams; tumor/normal paired
#### output: tumor/normal vcf files

### Lancet
[https://github.com/nygenome/lancet]  
Lancet uses a localized micro-assembly strategy to detect somatic mutation with high sensitivity and accuracy on a tumor/normal pair. Lancet is based on the colored de Bruijn graph assembly paradigm where tumor and normal reads are jointly analyzed within the same graph.

### Manta
[https://github.com/Illumina/manta]  
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a single efficient workflow.

### Mutect2
[https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2]  
[https://github.com/broadinstitute/gatk/releases]  
Call somatic short mutations via local assembly of haplotypes. Short mutations include SNVs and INDELs. 

### Strelka2
[https://github.com/Illumina/strelka]  
Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. The somatic calling model improves on the original Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell contamination in the normal sample. A final empirical variant re-scoring step using random forest models trained on various call quality features has been added to both callers to further improve precision.

### Varscan2
[https://github.com/dkoboldt/varscan]  
VarScan employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance.


### VarDict
[https://github.com/AstraZeneca-NGS/VarDictJava]  
VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability than many Java based variant callers.

***

## CopyNumber Variation

### Sequenza
[https://cran.r-project.org/web/packages/sequenza/index.html]  
Tools to analyze genomic sequencing data from paired normal-tumor samples, including cellularity and ploidy estimation; mutation and copy number (allele-specific and total copy number) detection, quantification and visualization.

### FACETS
[https://github.com/mskcc/facets]  
Algorithm to implement Fraction and Allele specific Copy number Estimate from Tumor/normal Sequencing.

***

## Somatic Variant Filtering
#### input: vcf file; one per variant caller
#### output: bcf file with "event" posterior probability values

### bcftools
  [https://github.com/samtools/bcftools]
  
  * bcftools concat > sort
    + Combine somatic snp/indel vcf data for each variant caller; position sort
  * bcftools isec
    + Produces an intersection matrix of the variants. This can be used to evaluate caller performance/metrics

### Varlociraptor
[https://varlociraptor.github.io/landing/]  
  * Calls SNVs, MNVs, indels, arbitrary replacements, inversions, duplications, haplotype blocks (combinations of any of the previous), and breakends.  
  * Supports all length ranges (from small to structural) with a unified statistical model.  
  * The statistical model entails all possible sources of uncertainty (mapping, typing, heterogeneity) and biases (strand, read pair orientation, read position, sampling, contamination, homologous regions).  
  * Resulting variant calls can be filtered by false discovery rate. No parameter tuning necessary.  
  * Maximum a posteriori allele frequency estimates are provided with each call.  

##### events
  1. germline
      + 'normal:0.5 | normal:1.0'
  2. somatic_tumor
      + 'normal:]0.0,0.5['
  3. somatic_normal
      + 'normal:0.0 & tumor:]0.0,1.0] & !$ffpe'
  * (optional) ffpe_artifact
      + '(C>T | G>A) & ((tumor:0.0 & normal:]0.0,0.1]) | (tumor:]0.0,0.1] & normal:0.0))'

***

## Annotation

### VEP

  * HGVS
  * ENSEMBL
  * dbSNP
  * SpliceAI
  * GNomAD
  * REVEL
  * ClinVAR
  * ~~LoFTee~~
  
## Additional Resources

### dbSNP
[https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/]  
Filtered for common snps.  
`bcftools filter -i'TYPE="snp" & COMMON=1'`  
