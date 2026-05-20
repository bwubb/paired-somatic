import csv

#ASCAT proto-pipeline

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

#with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
#ASCAT has dropped several samples.
#This probably wont be an issue in the future and should be changed back.
#We included pairs with germline >30X
with open('pair.table','r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

#I need FULL PATHS to bams
def bam_input(wildcards):
    bams=[]
    for t,n in PAIRS.items():
        bams+=[BAMS[t],BAMS[n]]
    return list(set(bams))

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

localrules:ASCAT_write_worksheet

#need segments file too for plotting purposes.
rule ascat_all:
    input:
        expand("data/work/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.report.csv",tumor=PAIRS.keys()),
        #expand("data/work/{tumor}/ascat/{tumor}.ascat.hrd.txt",tumor=PAIRS.keys())

rule ASCAT_write_worksheet:
    input:
        bam_input
    output:
        "data/work/ASCAT/worksheet.tsv"
    run:
        with open(output[0],'w') as f:
            writer=csv.writer(f,delimiter="\t")
            writer.writerow(['Patient_ID','Normal_ID','Normal_file','Gender'])
            for t,n in PAIRS.items():
                writer.writerow([t,n,BAMS[n],"XX"])


#output not directly used, picked X out of all any chr number.
rule ASCAT_prepareTargetedSeq:
    input:
        worksheet="data/work/ASCAT/worksheet.tsv"
    output:
        "data/work/ASCAT/alleleData/Cleaned/alleleData_chrX.txt",
        "data/work/ASCAT/alleleData/Cleaned/loci_chrX.txt"
    params:
        outdir="data/work/ASCAT",
        alleles=config['resources']['ascat_alleles'],#"G1000_allelesAll_hg38/G1000_alleles_hg38_chr"
        bed=config['resources']['targets_bed'],
        allelecounter="/home/bwubb/software/alleleCount/bin/alleleCounter",
        min_counts=10,
        min_base_qual=20,
        min_map_qual=35
    threads:
        8
    shell:
        "Rscript ASCAT_prepareTargetedSeq.R --worksheet {input.worksheet} --outdir {params.outdir} --alleles-prefix {params.alleles} --bed {params.bed} --allelecounter {params.allelecounter} --threads {threads} --min-counts {params.min_counts} --min-base-qual {params.min_base_qual} --min-map-qual {params.min_map_qual}"


rule ASCAT_prepareHTS:
    input:
        unpack(paired_bams),
        allele="data/work/ASCAT/alleleData/Cleaned/alleleData_chrX.txt",
        loci="data/work/ASCAT/alleleData/Cleaned/loci_chrX.txt"
    output:
        tumorLogR="data/work/{tumor}/ascat/Tumor_LogR.txt",
        tumorBAF="data/work/{tumor}/ascat/Tumor_BAF.txt",
        germlineLogR="data/work/{tumor}/ascat/Germline_LogR.txt",
        germlineBAF="data/work/{tumor}/ascat/Germline_BAF.txt"
    params:
        outdir="data/work/{tumor}/ascat",
        bed=config["resources"]["targets_bed"],
        alleles="data/work/ASCAT/alleleData/Cleaned/alleleData_chr",
        loci="data/work/ASCAT/alleleData/Cleaned/loci_chr",
        allelecounter="/home/bwubb/software/alleleCount/bin/alleleCounter",#need to fix the hardcoded path
        tumor='{tumor}',
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    threads: 1
    shell:
        "Rscript ASCAT_prepareHTS.R --tumor-bam {input.tumor} --normal-bam {input.normal} --tumor-name {params.tumor} --normal-name {params.normal} --outdir {params.outdir} --bed {params.bed} --alleles-prefix {params.alleles} --loci-prefix {params.loci} --allelecounter {params.allelecounter} --threads {threads}"

rule ASCAT_runAscat:
    input:
        tumorLogR="data/work/{tumor}/ascat/Tumor_LogR.txt",
        tumorBAF="data/work/{tumor}/ascat/Tumor_BAF.txt",
        germlineLogR="data/work/{tumor}/ascat/Germline_LogR.txt",
        germlineBAF="data/work/{tumor}/ascat/Germline_BAF.txt"
    output:
        segments="data/work/{tumor}/ascat/{tumor}.segments.txt",
        rdata="data/work/{tumor}/ascat/ASCAT_objects.Rdata"
    params:
        outdir="data/work/{tumor}/ascat",
        gc_file=config["resources"]["ascat_gc"]#gc_file="GC_G1000_hg38.txt"
    shell:
        "Rscript ASCAT_runAscat.R --tumor-logr {input.tumorLogR} --tumor-baf {input.tumorBAF} --germline-logr {input.germlineLogR} --germline-baf {input.germlineBAF} --outdir {params.outdir} --gc-file {params.gc_file} --rdata {output.rdata}"

rule ASCAT_to_bed:
    input:
        "data/work/{tumor}/ascat/{tumor}.segments.txt"
    output:
        "data/work/{tumor}/ascat/{tumor}.segments.bed"
    shell:
        "python cnv_to_bed.py -c ascat {input}"

rule ASCAT_annotsv:
    input:
        "data/work/{tumor}/ascat/{tumor}.segments.bed"
    output:
        "data/work/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule ASCAT_annotsv_parser:
    input:
        "data/work/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.tsv"
    output:
        "data/work/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule ascat_HRDex:
    input:
        "data/work/{tumor}/ascat/{tumor}.segments.bed"
    output:
        "data/work/{tumor}/ascat/{tumor}.ascat.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """
